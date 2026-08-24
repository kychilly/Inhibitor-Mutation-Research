#!/usr/bin/env python3
"""
summarize_lost_hbonds.py

Reads plip_analysis_results.json and exports detailed summary tables of both
hydrogen bond and hydrophobic contact changes (counts, residues, lost interactions)
directly to results/plip_results/hbond_loss_summary.txt.
"""

import json
from pathlib import Path

JSON_PATH = Path("results/plip_analysis/plip_analysis_results.json")
OUTPUT_DIR = Path("results/plip_results")
OUTPUT_FILE = OUTPUT_DIR / "hbond_loss_summary.txt"

DRUG_GENERATIONS = {
    "erlotinib": "Gen 1",
    "gefitinib": "Gen 1",
    "afatinib": "Gen 2",
    "dacomitinib": "Gen 2",
    "osimertinib": "Gen 3",
}


def build_hbond_summary(data):
    lines = []
    lines.append("=== HYDROGEN BOND LOSS SUMMARY (WILDTYPE VS T790M) ===")
    lines.append("")

    header = f"{'GEN':<6} | {'DRUG':<12} | {'WT COUNT':<9} | {'MUT COUNT':<9} | {'WT H-BONDS (RESIDUE & DIST)':<32} | {'LOST H-BONDS IN T790M':<25}"
    lines.append(header)
    lines.append("-" * len(header))

    for drug, gen in DRUG_GENERATIONS.items():
        wt_hb = data.get(f"wildtype_{drug}", {}).get("h_bonds", [])
        mut_hb = data.get(f"t790m_{drug}", {}).get("h_bonds", [])

        wt_count = len(wt_hb)
        mut_count = len(mut_hb)

        wt_str = ", ".join([f"{h['residue']} ({h['distance_angstrom']:.2f}Å)" for h in wt_hb]) if wt_hb else "None"

        mut_residues = {h["residue"] for h in mut_hb}
        lost_residues = [h["residue"] for h in wt_hb if h["residue"] not in mut_residues]
        lost_str = ", ".join(lost_residues) if lost_residues else (
            "None" if wt_count == mut_count else "Shifted Binding")

        line = f"{gen:<6} | {drug.capitalize():<12} | {wt_count:<9} | {mut_count:<9} | {wt_str:<32} | {lost_str:<25}"
        lines.append(line)

    return "\n".join(lines)


def build_hydrophobic_summary(data):
    lines = []
    lines.append("\n" + "=" * 110)
    lines.append("=== HYDROPHOBIC CONTACT LOSS SUMMARY (WILDTYPE VS T790M) ===")
    lines.append("Note: Hydrophobic contact tracking is not strictly required for primary H-bond loss objectives,")
    lines.append("      but is included here to provide a complete picture of non-covalent active site changes.")
    lines.append("=" * 110 + "\n")

    header = f"{'GEN':<6} | {'DRUG':<12} | {'WT COUNT':<9} | {'MUT COUNT':<9} | {'WT HYDROPHOBIC RESIDUES':<35} | {'LOST HYDROPHOBIC CONTACTS':<25}"
    lines.append(header)
    lines.append("-" * len(header))

    for drug, gen in DRUG_GENERATIONS.items():
        wt_hp = data.get(f"wildtype_{drug}", {}).get("hydrophobic", [])
        mut_hp = data.get(f"t790m_{drug}", {}).get("hydrophobic", [])

        wt_count = len(wt_hp)
        mut_count = len(mut_hp)

        wt_res_list = sorted(list({h["residue"] for h in wt_hp}))
        wt_str = ", ".join(wt_res_list) if wt_res_list else "None"

        mut_residues = {h["residue"] for h in mut_hp}
        lost_residues = [res for res in wt_res_list if res not in mut_residues]
        lost_str = ", ".join(lost_residues) if lost_residues else "None"

        line = f"{gen:<6} | {drug.capitalize():<12} | {wt_count:<9} | {mut_count:<9} | {wt_str:<35} | {lost_str:<25}"
        lines.append(line)

    return "\n".join(lines)


def main():
    if not JSON_PATH.exists():
        print(f"[-] Error: Could not find {JSON_PATH}")
        return

    with open(JSON_PATH, "r", encoding="utf-8") as f:
        data = json.load(f)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    hbond_section = build_hbond_summary(data)
    hydrophobic_section = build_hydrophobic_summary(data)

    full_output = f"{hbond_section}\n{hydrophobic_section}\n"

    with open(OUTPUT_FILE, "w", encoding="utf-8") as out_f:
        out_f.write(full_output)

    print(full_output)
    print(f"[+] Successfully saved combined summary to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()