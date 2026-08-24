#!/usr/bin/env python3
"""
verify_plip_results.py

Performs Step 2 (Active Site Residue Audit) and Step 3 (Spatial Geometry
Verification) directly on PLIP output files without requiring PyMOL.
"""

import json
import math
from pathlib import Path

JSON_PATH = Path("results/plip_analysis/plip_analysis_results.json")
PLIP_DIR = Path("results/plip_results")

# Key EGFR kinase domain residues expected in binding site interactions
KEY_EGFR_RESIDUES = {
    "MET793": "Hinge Region Anchor",
    "THR790": "Gatekeeper (WT)",
    "MET790": "Gatekeeper (T790M)",
    "CYS797": "Covalent Attachment Residue",
    "LYS745": "Catalytic Lysine",
    "ASP855": "DFG Motif Aspartate",
    "LEU718": "P-Loop",
    "PHE723": "P-Loop",
    "VAL726": "P-Loop",
    "ALA743": "Beta-sheet Core",
    "LEU844": "C-lobe Hydrophobic Pocket",
}


def verify_step_2_residues(data):
    """Step 2: Check for presence of key active site residues."""
    print("=== STEP 2: Active Site Residue Audit ===")

    for condition, details in data.items():
        found_key_residues = []

        all_interactions = details.get("h_bonds", []) + details.get("hydrophobic", [])
        for inter in all_interactions:
            res = inter["residue"]
            if res in KEY_EGFR_RESIDUES:
                found_key_residues.append(f"{res} ({KEY_EGFR_RESIDUES[res]})")

        found_key_residues = sorted(list(set(found_key_residues)))

        print(f"\n[+] Condition: {condition.upper()}")
        print(f"    Total Interactions: {details.get('total_interactions', 0)}")
        if found_key_residues:
            print("    Key EGFR Pocket Residues Hit:")
            for k_res in found_key_residues:
                print(f"      • {k_res}")
        else:
            print("    [!] Warning: No primary key pocket residues detected.")


def calculate_min_distance(complex_pdb):
    """Step 3: Calculate minimum Euclidean distance between Ligand and Receptor."""
    lig_coords = []
    rec_coords = []

    with open(complex_pdb, "w", encoding="utf-8") if not complex_pdb.exists() else open(complex_pdb, "r",
                                                                                        encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue

                resname = line[17:20].strip()
                if resname == "LIG":
                    lig_coords.append((x, y, z))
                else:
                    rec_coords.append((x, y, z))

    if not lig_coords or not rec_coords:
        return None

    min_sq_dist = float("inf")
    for lx, ly, lz in lig_coords:
        for rx, ry, rz in rec_coords:
            sq_dist = (lx - rx) ** 2 + (ly - ry) ** 2 + (lz - rz) ** 2
            if sq_dist < min_sq_dist:
                min_sq_dist = sq_dist

    return math.sqrt(min_sq_dist)


def verify_step_3_geometry(receptors, ligands):
    """Step 3: Geometry check across all generated PDB complexes."""
    print("\n\n=== STEP 3: Complex Geometry & Binding Distance Check ===")

    for rec in receptors:
        for lig in ligands:
            condition = f"{rec}_{lig}"
            complex_pdb = PLIP_DIR / rec / lig / f"{condition}_complex.pdb"

            if not complex_pdb.exists():
                # Check un-organized directory fallback
                complex_pdb = PLIP_DIR / f"{condition}_complex.pdb"

            if complex_pdb.exists():
                min_dist = calculate_min_distance(complex_pdb)
                if min_dist is not None:
                    status = "VALID (In Pocket)" if min_dist < 4.5 else "INVALID (Unbound/Floating)"
                    print(f"  • {condition:<20} | Min Dist: {min_dist:.2f} Å | Status: {status}")
                else:
                    print(f"  • {condition:<20} | Could not parse ligand/receptor coordinates")
            else:
                print(f"  • {condition:<20} | Complex PDB not found at {complex_pdb}")


def main():
    if not JSON_PATH.exists():
        print(f"Error: {JSON_PATH} not found.")
        return

    with open(JSON_PATH, "r", encoding="utf-8") as f:
        data = json.load(f)

    receptors = ["wildtype", "t790m"]
    ligands = ["afatinib", "dacomitinib", "erlotinib", "gefitinib", "osimertinib"]

    verify_step_2_residues(data)
    verify_step_3_geometry(receptors, ligands)


if __name__ == "__main__":
    main()