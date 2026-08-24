import os
import re
import json
import numpy as np
import pandas as pd
from pathlib import Path

# Base configuration
BASE_DIR = Path("data/processed/docking_results")

# Results config, this could really be a single line but oh well
Results_DIR = Path("results")

# Metadata mapping and chronological order
GENERATIONS = {
    "gefitinib": "1st Gen",
    "erlotinib": "1st Gen",
    "afatinib": "2nd Gen",
    "dacomitinib": "2nd Gen",
    "osimertinib": "3rd Gen"
}

ORDER_MAP = {
    "1st Gen": 1,
    "2nd Gen": 2,
    "3rd Gen": 3
}

SCORE_PATTERN = re.compile(r"^\s*1\s+(-?\d+\.\d+)")


def parse_raw_logs():
    """Reads all .log files and saves raw per-seed scores to CSV."""
    records = []
    for log_file in BASE_DIR.rglob("*.log"):
        target = log_file.parents[2].name.lower()
        ligand = log_file.parents[1].name.lower()

        seed_match = re.search(r"seed(\d+)", log_file.name)
        seed = int(seed_match.group(1)) if seed_match else None

        with open(log_file, "r") as f:
            for line in f:
                match = SCORE_PATTERN.match(line)
                if match:
                    records.append({
                        "target": target,
                        "ligand": ligand,
                        "seed": seed,
                        "affinity_kcal_mol": float(match.group(1)),
                        "log_file": str(log_file)
                    })
                    break

    df_raw = pd.DataFrame(records)
    out_csv = BASE_DIR / "raw_docking_scores.csv"
    df_raw.to_csv(out_csv, index=False)
    print(f"[+] Processed {len(df_raw)} log files -> {out_csv}")
    return df_raw


def generate_evaluation_json(df_raw):
    """Calculates grouped metrics, saves docking_evaluation.json, and returns chronological summary rows."""
    evaluation_results = {}
    summary_rows = []

    ligands = df_raw["ligand"].unique()

    for lig in sorted(ligands):
        lig_df = df_raw[df_raw["ligand"] == lig]
        gen = GENERATIONS.get(lig, "Unknown")

        wt_scores = lig_df[lig_df["target"] == "wildtype"]["affinity_kcal_mol"].tolist()
        t790m_scores = lig_df[lig_df["target"] == "t790m"]["affinity_kcal_mol"].tolist()

        wt_mean = float(np.mean(wt_scores)) if wt_scores else None
        wt_std = float(np.std(wt_scores, ddof=1)) if len(wt_scores) > 1 else 0.0

        t790m_mean = float(np.mean(t790m_scores)) if t790m_scores else None
        t790m_std = float(np.std(t790m_scores, ddof=1)) if len(t790m_scores) > 1 else 0.0

        ddg = (t790m_mean - wt_mean) if (t790m_mean is not None and wt_mean is not None) else None

        # Build JSON data structure
        evaluation_results[lig] = {
            "generation": gen,
            "wildtype": {
                "raw_scores": wt_scores,
                "mean": round(wt_mean, 3) if wt_mean is not None else None,
                "std": round(wt_std, 3),
                "formatted": f"{wt_mean:.2f} ± {wt_std:.2f}" if wt_mean is not None else "N/A"
            },
            "t790m": {
                "raw_scores": t790m_scores,
                "mean": round(t790m_mean, 3) if t790m_mean is not None else None,
                "std": round(t790m_std, 3),
                "formatted": f"{t790m_mean:.2f} ± {t790m_std:.2f}" if t790m_mean is not None else "N/A"
            },
            "ddG_shift": {
                "val": round(ddg, 3) if ddg is not None else None,
                "formatted": f"{ddg:+.2f}" if ddg is not None else "N/A"
            }
        }

        # Build tabular row object
        summary_rows.append({
            "Ligand": lig.capitalize(),
            "Target": gen,
            "WT Mean ± SD": f"{wt_mean:.2f} ± {wt_std:.2f}" if wt_mean is not None else "N/A",
            "T790M Mean ± SD": f"{t790m_mean:.2f} ± {t790m_std:.2f}" if t790m_mean is not None else "N/A",
            "Score Shift (ΔΔG)": f"{ddg:+.2f}" if ddg is not None else "N/A",
            "sort_key": ORDER_MAP.get(gen, 99)
        })

    # Save JSON file
    json_path = BASE_DIR / "docking_evaluation.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(evaluation_results, f, indent=4)
    print(f"[+] Saved evaluation metric JSON -> {json_path}")

    # Sort chronologically (1st Gen -> 2nd Gen -> 3rd Gen), then alphabetically
    sorted_rows = sorted(summary_rows, key=lambda x: (x["sort_key"], x["Ligand"]))
    for r in sorted_rows:
        del r["sort_key"]

    return sorted_rows


def export_summary_tables(summary_rows):
    """Exports structured tables in Markdown, CSV, and LaTeX formats."""
    df_summary = pd.DataFrame(summary_rows)

    # 1. Export CSV Summary
    csv_path = BASE_DIR / "summary_table.csv"
    df_summary.to_csv(csv_path, index=False)

    # 2. Export Markdown Summary
    md_path = BASE_DIR / "summary_table.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(df_summary.to_markdown(index=False))

    # 3. Export LaTeX Summary
    tex_path = BASE_DIR / "summary_table.tex"
    latex_content = (
        "\\begin{table}[htbp]\n"
        "\\centering\n"
        "\\caption{Binding Affinities (kcal/mol) and Mutation Score Shifts Across EGFR Generations}\n"
        "\\label{tab:docking_summary}\n"
        "\\begin{tabular}{lcccc}\n"
        "\\toprule\n"
        "Ligand & Generation & WT Mean $\\pm$ SD & T790M Mean $\\pm$ SD & Score Shift ($\\Delta\\Delta G$) \\\\\n"
        "\\midrule\n"
    )
    for r in summary_rows:
        latex_content += f"{r['Ligand']} & {r['Target']} & {r['WT Mean ± SD']} & {r['T790M Mean ± SD']} & {r['Score Shift (ΔΔG)']} \\\\\n"

    latex_content += (
        "\\bottomrule\n"
        "\\end{tabular}\n"
        "\\end{table}\n"
    )

    with open(tex_path, "w", encoding="utf-8") as f:
        f.write(latex_content)

    print(f"[+] Saved Summary CSV    -> {csv_path}")
    print(f"[+] Saved Summary MD     -> {md_path}")
    print(f"[+] Saved Summary LaTeX  -> {tex_path}")

    # Just storing this same result in the results folder as well yk
    Path("results").mkdir(exist_ok=True)
    df_summary.to_csv("results/docking_table.csv", index=False)


if __name__ == "__main__":
    df_raw = parse_raw_logs()
    summary_rows = generate_evaluation_json(df_raw)
    export_summary_tables(summary_rows)