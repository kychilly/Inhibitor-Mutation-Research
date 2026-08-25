#!/usr/bin/env python3
"""
build_results_table.py

Aggregates Vina docking statistics (Top ΔG, Mean ΔG, Std Dev) and PLIP non-covalent
interaction counts across all 10 conditions. Exports the final summary table directly
to Markdown and LaTeX formats in results/final_summary_table.
"""

import json
import numpy as np
from pathlib import Path

# File paths
RESULTS_DIR = Path("data/processed/docking_results")
PLIP_JSON = Path("results/plip_analysis/plip_analysis_results.json")
OUTPUT_DIR = Path("results/final_summary_table")

DRUG_GENERATIONS = {
    "erlotinib": "Gen 1",
    "gefitinib": "Gen 1",
    "afatinib": "Gen 2",
    "dacomitinib": "Gen 2",
    "osimertinib": "Gen 3",
}

RECEPTORS = {
    "wildtype": "Wildtype (WT)",
    "t790m": "T790M Mutant",
}


def parse_vina_energies(rec, lig):
    """Scans all seed*.pdbqt files for a condition to calculate Top ΔG, Mean ΔG, and Std Dev."""
    pdbqt_dir = RESULTS_DIR / rec / lig / "pdbqt"
    candidates = sorted(pdbqt_dir.glob(f"{rec}_{lig}_seed*.pdbqt"))

    energies = []
    for path in candidates:
        try:
            with open(path, "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("REMARK VINA RESULT:"):
                        affinity = float(line.split()[3])
                        energies.append(affinity)
                        break
        except (OSError, IndexError, ValueError):
            continue

    if not energies:
        return None, None, None

    top_dg = min(energies)
    mean_dg = np.mean(energies)
    std_dg = np.std(energies)

    return top_dg, mean_dg, std_dg


def generate_markdown(rows):
    """Generates clean, column-aligned Markdown table with matching 9 headers and alignments."""
    md = []
    md.append("# Summary of Docking Affinities and PLIP Interaction Profiling\n")

    # 9 Headers
    header = "| Target | Variant       | Drug        | Generation | Top ΔG (kcal/mol) | Mean ΔG ± SD (kcal/mol) | H-Bonds | Hydrophobic | Total Interactions |"
    divider = "| :---   | :---          | :---        | :---:      | :---:             | :---:                   | :---:   | :---:       | :---:              |"

    md.append(header)
    md.append(divider)

    for r in rows:
        target_str = f"{r['target']:<6}"
        variant_str = f"{r['variant']:<13}"
        drug_str = f"{r['drug']:<11}"
        gen_str = f"{r['gen']:<10}"
        top_str = f"{r['top_dg']:-17.2f}"
        mean_str = f"{r['mean_dg']:.2f} ± {r['std_dg']:.2f}".center(23)
        hb_str = f"{r['h_bonds']:^7}"
        hp_str = f"{r['hydrophobic']:^11}"
        tot_str = f"{r['total_int']:^18}"

        md.append(
            f"| {target_str} | {variant_str} | {drug_str} | {gen_str} | {top_str} | {mean_str} | {hb_str} | {hp_str} | {tot_str} |")

    return "\n".join(md)


def generate_latex(rows):
    """Generates camera-ready LaTeX table for manuscript insertion."""
    tex = []
    tex.append("% LaTeX Table for Results Section")
    tex.append("\\begin{table}[htbp]")
    tex.append("  \\centering")
    tex.append(
        "  \\caption{Summary of AutoDock Vina binding energies and PLIP interaction counts across EGFR Wildtype and T790M variants.}")
    tex.append("  \\label{tab:docking_plip_summary}")
    tex.append("  \\begin{tabular}{lllcS[table-format=-1.2]cSSS}")
    tex.append("    \\toprule")
    tex.append(
        "    \\textbf{Target} & \\textbf{Variant} & \\textbf{Drug} & \\textbf{Gen.} & {\\textbf{Top $\\Delta$G (kcal/mol)}} & {\\textbf{Mean $\\Delta$G $\\pm$ SD}} & {\\textbf{H-Bonds}} & {\\textbf{Hydrophobic}} & {\\textbf{Total}} \\\\")
    tex.append("    \\midrule")

    current_variant = None
    for r in rows:
        if current_variant and current_variant != r['variant']:
            tex.append("    \\midrule")
        current_variant = r['variant']

        tex.append(
            f"    EGFR & {r['variant']} & {r['drug']} & {r['gen']} & {r['top_dg']:.2f} & "
            f"{r['mean_dg']:.2f} $\\pm$ {r['std_dg']:.2f} & {r['h_bonds']} & {r['hydrophobic']} & {r['total_int']} \\\\"
        )

    tex.append("    \\bottomrule")
    tex.append("  \\end{tabular}")
    tex.append("\\end{table}")

    return "\n".join(tex)


def main():
    if not PLIP_JSON.exists():
        print(f"[-] Error: Could not find PLIP JSON at {PLIP_JSON}")
        return

    with open(PLIP_JSON, "r", encoding="utf-8") as f:
        plip_data = json.load(f)

    rows = []

    for rec_key, rec_label in RECEPTORS.items():
        for lig_key, gen in DRUG_GENERATIONS.items():
            cond_key = f"{rec_key}_{lig_key}"

            top_dg, mean_dg, std_dg = parse_vina_energies(rec_key, lig_key)

            cond_plip = plip_data.get(cond_key, {})
            h_bonds = len(cond_plip.get("h_bonds", []))
            hydrophobic = len(cond_plip.get("hydrophobic", []))
            total_int = cond_plip.get("total_interactions", h_bonds + hydrophobic)

            rows.append({
                "target": "EGFR",
                "variant": rec_label,
                "drug": lig_key.capitalize(),
                "gen": gen,
                "top_dg": top_dg if top_dg is not None else 0.0,
                "mean_dg": mean_dg if mean_dg is not None else 0.0,
                "std_dg": std_dg if std_dg is not None else 0.0,
                "h_bonds": h_bonds,
                "hydrophobic": hydrophobic,
                "total_int": total_int,
            })

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    md_content = generate_markdown(rows)
    latex_content = generate_latex(rows)

    (OUTPUT_DIR / "summary_table.md").write_text(md_content, encoding="utf-8")
    (OUTPUT_DIR / "summary_table.tex").write_text(latex_content, encoding="utf-8")

    print(md_content)
    print(f"\n[+] Saved Markdown table to: {OUTPUT_DIR / 'summary_table.md'}")
    print(f"[+] Saved LaTeX table to:    {OUTPUT_DIR / 'summary_table.tex'}")


if __name__ == "__main__":
    main()