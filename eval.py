import re
import json
import numpy as np
import pandas as pd
from pathlib import Path

# Paths
BASE_DIR = Path("data/processed/docking_results")
JSON_OUT = BASE_DIR / "docking_evaluation.json"
MD_OUT = BASE_DIR / "docking_evaluation_report.md"

# Generation metadata mapping
GENERATIONS = {
    "gefitinib": "1st Gen",
    "erlotinib": "1st Gen",
    "afatinib": "2nd Gen",
    "dacomitinib": "2nd Gen",
    "osimertinib": "3rd Gen"
}

# Chronological sorting order
ORDER_MAP = {
    "1st Gen": 1,
    "2nd Gen": 2,
    "3rd Gen": 3
}

# Regex pattern for Vina mode 1 score
SCORE_PATTERN = re.compile(r"^\s*1\s+(-?\d+\.\d+)")


def parse_and_evaluate():
    raw_data = []

    # 1. Parse all log files
    for log_file in BASE_DIR.rglob("*.log"):
        target = log_file.parents[2].name.lower()
        ligand = log_file.parents[1].name.lower()

        seed_match = re.search(r"seed(\d+)", log_file.name)
        seed = int(seed_match.group(1)) if seed_match else None

        with open(log_file, "r") as f:
            for line in f:
                match = SCORE_PATTERN.match(line)
                if match:
                    raw_data.append({
                        "target": target,
                        "ligand": ligand,
                        "seed": seed,
                        "affinity": float(match.group(1))
                    })
                    break

    # 2. Convert to DataFrame for aggregation
    df = pd.DataFrame(raw_data)

    # 3. Calculate summary metrics per ligand & target
    evaluation_results = {}
    ligands = df["ligand"].unique()

    for lig in ligands:
        lig_df = df[df["ligand"] == lig]

        wt_scores = lig_df[lig_df["target"] == "wildtype"]["affinity"].tolist()
        t790m_scores = lig_df[lig_df["target"] == "t790m"]["affinity"].tolist()

        wt_mean = float(np.mean(wt_scores)) if wt_scores else None
        wt_std = float(np.std(wt_scores, ddof=1)) if len(wt_scores) > 1 else 0.0

        t790m_mean = float(np.mean(t790m_scores)) if t790m_scores else None
        t790m_std = float(np.std(t790m_scores, ddof=1)) if len(t790m_scores) > 1 else 0.0

        # Delta Delta G Shift = T790M Mean - WT Mean
        ddg = (t790m_mean - wt_mean) if (t790m_mean is not None and wt_mean is not None) else None

        evaluation_results[lig] = {
            "generation": GENERATIONS.get(lig, "Unknown"),
            "wildtype": {
                "raw_scores": wt_scores,
                "mean": round(wt_mean, 3) if wt_mean is not None else None,
                "std": round(wt_std, 3),
                "formatted": f"{wt_mean:.2f} ± {wt_std:.2f}" if wt_mean is not None else "N/A"
            },
            "t790m": {
                "raw_scores": t790m_scores,
                "mean": round(t790m_std, 3),
                "std": round(t790m_std, 3),
                "formatted": f"{t790m_mean:.2f} ± {t790m_std:.2f}" if t790m_mean is not None else "N/A"
            },
            "ddG_shift": {
                "val": round(ddg, 3) if ddg is not None else None,
                "formatted": f"{ddg:+.2f}" if ddg is not None else "N/A"
            }
        }

    # 4. Save to JSON
    with open(JSON_OUT, "w", encoding="utf-8") as f:
        json.dump(evaluation_results, f, indent=4)

    print(f"[+] Successfully evaluated docking scores and saved JSON to: {JSON_OUT}")
    return evaluation_results


def generate_markdown_report(json_path):
    """Loads JSON data and builds a chronologically sorted Markdown report file."""
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    rows = []
    for lig, info in data.items():
        gen = info.get("generation", "Unknown")
        rows.append({
            "Ligand": lig.capitalize(),
            "Target": gen,
            "WT Mean ± SD": info["wildtype"]["formatted"],
            "T790M Mean ± SD": info["t790m"]["formatted"],
            "Score Shift (ΔΔG)": info["ddG_shift"]["formatted"],
            "sort_key": ORDER_MAP.get(gen, 99)
        })

    # Sort chronologically by generation, then alphabetically by ligand name
    sorted_rows = sorted(rows, key=lambda x: (x["sort_key"], x["Ligand"]))
    for r in sorted_rows:
        del r["sort_key"]

    df = pd.DataFrame(sorted_rows)

    md_content = f"""# Docking Evaluation Summary Report

{df.to_markdown(index=False)}

---
* **Note on ΔΔG:** Calculated as `T790M Mean - WT Mean`. Positive values represent an energetic penalty against the mutation, whereas negative values indicate stronger relative computed affinity.
"""

    with open(MD_OUT, "w", encoding="utf-8") as f:
        f.write(md_content)

    print(f"[+] Saved human-readable Markdown evaluation to: {MD_OUT}\n")
    return sorted_rows


def print_summary_table(rows):
    """Printable terminal summary directly from sorted data."""
    headers = ["Ligand", "Target", "WT Mean ± SD", "T790M Mean ± SD", "Score Shift (ΔΔG)"]
    print(f"{headers[0]:<15} {headers[1]:<10} {headers[2]:<18} {headers[3]:<18} {headers[4]:<18}")
    print("-" * 82)

    for row in rows:
        print(f"{row['Ligand']:<15} {row['Target']:<10} {row['WT Mean ± SD']:<18} {row['T790M Mean ± SD']:<18} {row['Score Shift (ΔΔG)']:<18}")


if __name__ == "__main__":
    parse_and_evaluate()
    sorted_rows = generate_markdown_report(JSON_OUT)
    print_summary_table(sorted_rows)