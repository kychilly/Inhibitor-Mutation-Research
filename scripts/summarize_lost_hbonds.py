#!/usr/bin/env python3
"""
Parses PLIP XML results, tracks unique residue losses/gains,
and writes results to results/hbond_loss_summary.txt.
"""

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

PLIP_DIR = Path("results/plip_results")
OUTPUT_FILE = Path("results/hbond_loss_summary.txt")

DRUG_GENS = {
    "erlotinib": ("Gen 1", "Erlotinib"),
    "gefitinib": ("Gen 1", "Gefitinib"),
    "afatinib": ("Gen 2", "Afatinib"),
    "dacomitinib": ("Gen 2", "Dacomitinib"),
    "osimertinib": ("Gen 3", "Osimertinib"),
}


def find_xml_file(state_folder: str, drug_key: str) -> Path | None:
    target_dir = PLIP_DIR / state_folder / drug_key
    if not target_dir.exists():
        return None
    matches = list(target_dir.glob(f"{state_folder}_{drug_key}*.xml")) or list(
        target_dir.glob("*.xml")
    )
    return matches[0] if matches else None


def extract_distance(hb_elem: ET.Element) -> str:
    for tag in ["dist_h_a", "dist_d_a", "dist"]:
        val = hb_elem.findtext(tag)
        if val and val.strip() and val.strip() != "0.0":
            try:
                return f"{float(val.strip()):.2f}"
            except ValueError:
                return val.strip()
    return "N/A"


def parse_plip_xml(xml_path: Path):
    h_bonds = []
    hydrophobics = []

    if not xml_path or not xml_path.exists():
        return {"h_bonds": h_bonds, "hydrophobic": hydrophobics}

    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except ET.ParseError:
        return {"h_bonds": h_bonds, "hydrophobic": hydrophobics}

    for site in root.findall(".//bindingsite"):
        for hb in site.findall(".//hydrogen_bond"):
            resnr = hb.findtext("resnr", "").strip()
            restype = hb.findtext("restype", "").strip()
            dist_str = extract_distance(hb)
            res_label = f"{restype}{resnr}"
            h_bonds.append(
                {
                    "res": res_label,
                    "dist": dist_str,
                    "str": (
                        f"{res_label} ({dist_str}Å)"
                        if dist_str != "N/A"
                        else res_label
                    ),
                }
            )

        for hp in site.findall(".//hydrophobic_interaction"):
            resnr = hp.findtext("resnr", "").strip()
            restype = hp.findtext("restype", "").strip()
            hydrophobics.append({"res": f"{restype}{resnr}"})

    return {"h_bonds": h_bonds, "hydrophobic": hydrophobics}


def main():
    if not PLIP_DIR.exists():
        sys.exit(1)

    hb_rows = []
    hp_rows = []

    for drug_key, (gen, drug_name) in DRUG_GENS.items():
        wt_xml = find_xml_file("wildtype", drug_key)
        mut_xml = find_xml_file("t790m", drug_key)

        wt_data = parse_plip_xml(wt_xml)
        mut_data = parse_plip_xml(mut_xml)

        # H-Bonds processing
        wt_hb = wt_data["h_bonds"]
        mut_hb = mut_data["h_bonds"]

        wt_hb_res = list(dict.fromkeys([h["res"] for h in wt_hb]))
        mut_hb_res = list(dict.fromkeys([h["res"] for h in mut_hb]))

        wt_hb_str = ", ".join([h["str"] for h in wt_hb]) if wt_hb else "None"

        lost_hb = [res for res in wt_hb_res if res not in mut_hb_res]
        gained_hb = [res for res in mut_hb_res if res not in wt_hb_res]

        lost_str = ", ".join(lost_hb) if lost_hb else "None"
        if gained_hb:
            lost_str += f" (Gained: {', '.join(gained_hb)})"

        hb_rows.append(
            (gen, drug_name, len(wt_hb), len(mut_hb), wt_hb_str, lost_str)
        )

        # Hydrophobic processing
        wt_hp = wt_data["hydrophobic"]
        mut_hp = mut_data["hydrophobic"]

        wt_hp_res = sorted(list(set([h["res"] for h in wt_hp])))
        mut_hp_res = set([h["res"] for h in mut_hp])

        wt_hp_str = ", ".join(wt_hp_res) if wt_hp_res else "None"

        lost_hp = [res for res in wt_hp_res if res not in mut_hp_res]
        gained_hp = sorted(
            [res for res in mut_hp_res if res not in wt_hp_res]
        )

        lost_hp_str = ", ".join(lost_hp) if lost_hp else "None"
        if gained_hp:
            lost_hp_str += f" (Gained: {', '.join(gained_hp)})"

        hp_rows.append(
            (gen, drug_name, len(wt_hp), len(mut_hp), wt_hp_str, lost_hp_str)
        )

    out = []
    out.append("=== HYDROGEN BOND LOSS SUMMARY (WILDTYPE VS T790M) ===\n")
    out.append(
        f"{'GEN':<6} | {'DRUG':<12} | {'WT COUNT':<9} | {'MUT COUNT':<9} | {'WT H-BONDS (RESIDUE & DIST)':<36} | {'LOST H-BONDS IN T790M (AND GAINED)':<40}"
    )
    out.append("-" * 128)
    for row in hb_rows:
        out.append(
            f"{row[0]:<6} | {row[1]:<12} | {row[2]:<9} | {row[3]:<9} | {row[4]:<36} | {row[5]:<40}"
        )

    out.append("\n\n" + "=" * 128)
    out.append("=== HYDROPHOBIC CONTACT LOSS SUMMARY (WILDTYPE VS T790M) ===")
    out.append("=" * 128 + "\n")
    out.append(
        f"{'GEN':<6} | {'DRUG':<12} | {'WT COUNT':<9} | {'MUT COUNT':<9} | {'WT HYDROPHOBIC RESIDUES':<45} | {'LOST / GAINED CONTACTS IN T790M':<45}"
    )
    out.append("-" * 139)
    for row in hp_rows:
        out.append(
            f"{row[0]:<6} | {row[1]:<12} | {row[2]:<9} | {row[3]:<9} | {row[4]:<45} | {row[5]:<45}"
        )

    final_text = "\n".join(out)
    print(final_text)

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text(final_text, encoding="utf-8")
    print(f"\n[+] Updated file: {OUTPUT_FILE.resolve()}")


if __name__ == "__main__":
    main()