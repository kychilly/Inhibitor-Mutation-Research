#!/usr/bin/env python3
"""
Automates PyMOL rendering for manuscript wildtype figures:
  - Figure 2C (or WT equivalent): Wildtype EGFR (T790) complexed with Erlotinib
               showing unhindered binding and ample spatial clearance.
  - Figure 2D (or WT equivalent): Wildtype EGFR (T790) complexed with Osimertinib
               showing baseline non-covalent active site orientation.
"""

import sys
from pathlib import Path

import numpy as np

try:
    import pymol
    from pymol import cmd
except ImportError:
    print("[-] Error: PyMOL module not found. Please run this inside an environment with PyMOL installed.")
    sys.exit(1)

OUTPUT_DIR = Path("results/figures")
PLIP_DIR = Path("results/plip_results/wildtype")  # Directory for WT complexes

# Distinct, high-contrast palette matching the mutant script pipeline
COLOR_CARTOON = "skyblue"       # protein backbone - clearly blue, not gray
COLOR_POCKET = "yellow"         # nearby side chains for context - clearly yellow
COLOR_T790_A = "deepteal"       # wildtype Thr790 residue for Erlotinib panel
COLOR_T790_B = "forest"         # wildtype Thr790 residue for Osimertinib panel
COLOR_LIGAND_A = "cyan"         # Erlotinib - bright cyan
COLOR_LIGAND_B = "orange"       # Osimertinib - bright orange


def find_pdb_file(drug_name):
    drug_dir = PLIP_DIR / drug_name

    # 1. Look inside the drug's dedicated directory (results/plip_results/wildtype/<drug_name>/)
    if drug_dir.exists():
        # Match "wildtype" or "wt" prefix dynamically
        candidates = sorted(list(drug_dir.glob(f"plipfixed.*{drug_name}_complex_*.pdb")))
        if not candidates:
            # Fallback to standard complex PDB if plipfixed is missing
            candidates = sorted(list(drug_dir.glob(f"*{drug_name}_complex*.pdb")))
    else:
        # 2. Search root wildtype folder if subdirectory does not exist
        candidates = sorted(list(PLIP_DIR.glob(f"plipfixed.*{drug_name}_complex_*.pdb")))
        if not candidates:
            candidates = sorted(list(PLIP_DIR.glob(f"*{drug_name}_complex*.pdb")))

    if not candidates:
        raise FileNotFoundError(f"Could not locate Wildtype PDB file for '{drug_name}' in {PLIP_DIR}")

    # Explicitly print which file is chosen so you can verify it match your desired model
    chosen = candidates[0]
    print(f"  [+] Found PDB for {drug_name}: {chosen.relative_to(PLIP_DIR)}")
    return chosen


def closest_atom_pair(sel1, sel2):
    """Finds the single closest heavy-atom pair between two selections and
    returns (distance_in_angstroms, atom1, atom2)."""
    model1 = cmd.get_model(sel1)
    model2 = cmd.get_model(sel2)
    if not model1.atom or not model2.atom:
        return None

    coords1 = np.array([a.coord for a in model1.atom])
    coords2 = np.array([a.coord for a in model2.atom])
    diffs = coords1[:, None, :] - coords2[None, :, :]
    dists = np.linalg.norm(diffs, axis=2)
    i, j = np.unravel_index(np.argmin(dists), dists.shape)
    return float(dists[i, j]), model1.atom[i], model2.atom[j]


def setup_pymol_display():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("stick_radius", 0.20)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.set("surface_quality", 2)
    cmd.set("transparency_mode", 2)
    cmd.set("depth_cue", 0)


def render_wt_pocket_figure(pdb_file, output_png, ligand_resn_list, t790_color, ligand_color, label,
                             label_offset_mode="camera_up"):
    """Shared rendering logic for Wildtype figures targeting Thr790."""
    print(f"[+] Rendering {label} using {pdb_file.name}...")
    setup_pymol_display()

    cmd.load(str(pdb_file), "complex")
    obj = "complex"

    resn_selector = " or ".join(f"resn {r}" for r in ligand_resn_list)
    cmd.select("ligand", f"{obj} and ({resn_selector})")
    cmd.select("t790", f"{obj} and resi 790")

    if cmd.count_atoms("ligand") == 0:
        print(f"  [!] Warning: no ligand atoms matched {ligand_resn_list} - check residue naming in the PDB.")
    if cmd.count_atoms("t790") == 0:
        print("  [!] Warning: no atoms matched resi 790 - check residue numbering in the PDB.")

    # Measure the actual closest heavy-atom distance between T790 (Threonine) and ligand
    result = closest_atom_pair("t790", "ligand")
    if result is not None:
        min_dist, atom_t790, atom_ligand = result
        print(f"  [i] Closest heavy-atom contact: T790/{atom_t790.name} <-> "
              f"{atom_ligand.resn}/{atom_ligand.name}  =  {min_dist:.2f} A")

        cmd.select("closest_t790", f"t790 and id {atom_t790.id}")
        cmd.select("closest_ligand", f"ligand and id {atom_ligand.id}")
        cmd.select("contact_t790", "t790 within 4.5 of ligand")
        cmd.select("contact_ligand", "ligand within 4.5 of t790")
        contact_atoms = cmd.count_atoms("contact_t790") + cmd.count_atoms("contact_ligand")
    else:
        min_dist = None
        contact_atoms = 0
        print("  [!] Could not measure distance - t790 or ligand selection is empty.")

    # Background context residues
    cmd.select("pocket", f"byres ({obj} within 5.0 of (ligand or t790)) and not t790")

    cmd.hide("everything")

    # Cartoon: transparent blue
    cmd.show("cartoon", obj)
    cmd.color(COLOR_CARTOON, obj)
    cmd.set("cartoon_transparency", 0.75, obj)

    # Pocket context residues: thin yellow sticks
    cmd.show("sticks", "pocket")
    cmd.color(COLOR_POCKET, "pocket and elem C")
    cmd.set("stick_radius", 0.10, "pocket")

    # T790 (Threonine): solid color + surface representation
    cmd.show("sticks", "t790")
    cmd.color(t790_color, "t790 and elem C")
    cmd.show("surface", "t790")
    cmd.color(t790_color, "t790")
    cmd.set("transparency", 0.2, "t790")

    # Ligand: solid color + surface representation
    cmd.show("sticks", "ligand")
    cmd.color(ligand_color, "ligand and elem C")
    cmd.set("stick_radius", 0.24, "ligand")
    cmd.show("surface", "ligand")
    cmd.color(ligand_color, "ligand")
    cmd.set("transparency", 0.45, "ligand")

    # Draw measured distance line
    if min_dist is not None:
        cmd.distance("contact_dist", "closest_t790", "closest_ligand")
        cmd.hide("labels", "contact_dist")
        cmd.color("black", "contact_dist")
        cmd.set("dash_width", 4, "contact_dist")
        cmd.set("dash_gap", 0.3, "contact_dist")
        cmd.set("dash_radius", 0.04, "contact_dist")

    # Camera framing
    if min_dist is not None:
        cmd.orient("closest_t790 or closest_ligand")
        cmd.zoom("closest_t790 or closest_ligand", buffer=6.0)
    elif contact_atoms > 0:
        cmd.orient("contact_t790 or contact_ligand")
        cmd.zoom("contact_t790 or contact_ligand", buffer=3.0)
    else:
        cmd.orient("t790 or ligand")
        cmd.zoom("t790 or ligand", buffer=4.0)

    cmd.turn("y", 15)

    # Label placement
    if min_dist is not None:
        ligand_coord = np.array(cmd.get_atom_coords("closest_ligand"))

        if label_offset_mode == "fixed":
            anchor_pos = ligand_coord + np.array([0.0, 2.5, 0.0])
        else:  # "camera_up"
            view = cmd.get_view()
            rot = view[0:9]
            up = np.array([rot[1], rot[4], rot[7]])
            norm = np.linalg.norm(up)
            if norm > 0:
                up = up / norm
            anchor_pos = ligand_coord + up * 3.5

        cmd.pseudoatom("dist_label_anchor", pos=list(anchor_pos))
        cmd.hide("everything", "dist_label_anchor")
        cmd.set("label_size", 32)
        cmd.set("label_color", "black")
        cmd.set("label_outline_color", "white")
        cmd.set("label_font_id", 7)
        cmd.label("dist_label_anchor", f'"{min_dist:.2f} A"')
        cmd.show("labels", "dist_label_anchor")

    cmd.ray(2400, 1800)
    cmd.png(str(output_png), dpi=300)
    print(f"  [Saved] {label} -> {output_png}")


def render_figure_wt_erlotinib(erlotinib_pdb, output_png):
    render_wt_pocket_figure(
        erlotinib_pdb, output_png,
        ligand_resn_list=["UNK", "LIG", "ERL"],
        t790_color=COLOR_T790_A,
        ligand_color=COLOR_LIGAND_A,
        label="Figure WT-Erlotinib (WT baseline clearance)",
        label_offset_mode="fixed",
    )


def render_figure_wt_osimertinib(osimertinib_pdb, output_png):
    render_wt_pocket_figure(
        osimertinib_pdb, output_png,
        ligand_resn_list=["UNK", "LIG", "OSI"],
        t790_color=COLOR_T790_B,
        ligand_color=COLOR_LIGAND_B,
        label="Figure WT-Osimertinib (WT baseline clearance)",
        label_offset_mode="camera_up",
    )


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    try:
        erlotinib_pdb = find_pdb_file("erlotinib")
        osimertinib_pdb = find_pdb_file("osimertinib")
    except FileNotFoundError as err:
        print(f"[-] File Resolution Error: {err}")
        sys.exit(1)

    wt_erl_path = OUTPUT_DIR / "Figure_WT_Erlotinib_Clearance.png"
    wt_osi_path = OUTPUT_DIR / "Figure_WT_Osimertinib_Clearance.png"

    render_figure_wt_erlotinib(erlotinib_pdb, wt_erl_path)
    render_figure_wt_osimertinib(osimertinib_pdb, wt_osi_path)

    print("\n[+] Both Wildtype PyMOL figures successfully rendered and exported to results/figures/")


if __name__ == "__main__":
    main()