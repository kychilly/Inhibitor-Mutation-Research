#!/usr/bin/env python3
"""
Automates PyMOL rendering for manuscript figures:
  - Figure 2A: T790M mutant complexed with Erlotinib showing steric collision
               between the bulky Methionine 790 side-chain and Erlotinib's aniline ring.
  - Figure 2B: T790M mutant complexed with Osimertinib showing structural
               accommodation of the M790 side chain by Osimertinib's flexible scaffold.
"""

# For figure 2A, just use scaryThing.png, since it puts the affinity number more clean
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
PLIP_DIR = Path("results/plip_results/t790m")

# Distinct, high-contrast palette so nothing reads as "gray blending into gray"
COLOR_CARTOON = "skyblue"       # protein backbone - clearly blue, not gray
COLOR_POCKET = "yellow"         # nearby side chains for context - clearly yellow
COLOR_M790_A = "red"            # mutated residue in fig 2A - solid red
COLOR_M790_B = "green"          # mutated residue in fig 2B - solid green
COLOR_LIGAND_A = "cyan"         # Erlotinib - bright cyan (distinct from red/yellow/blue)
COLOR_LIGAND_B = "orange"       # Osimertinib - bright orange (distinct from green/yellow/blue)


def find_pdb_file(drug_name):
    drug_dir = PLIP_DIR / drug_name
    if not drug_dir.exists():
        candidates = list(PLIP_DIR.glob(f"plipfixed.t790m_{drug_name}_complex_*.pdb"))
    else:
        candidates = list(drug_dir.glob(f"plipfixed.t790m_{drug_name}_complex_*.pdb"))

    if not candidates:
        raise FileNotFoundError(f"Could not locate PDB file for T790M + {drug_name} in {PLIP_DIR}")

    return candidates[0]


def closest_atom_pair(sel1, sel2):
    """Finds the single closest heavy-atom pair between two selections and
    returns (distance_in_angstroms, atom1, atom2). This is the actual
    measurement that matters for a 'steric clash' claim - '5 A between atom
    centers' is NOT the same as atoms touching, since each atom has its own
    ~1.5-1.8 A van der Waals radius. True clashes/close contacts are
    typically in the ~3.0-3.8 A range between heavy atom centers."""
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


def render_pocket_figure(pdb_file, output_png, ligand_resn_list, m790_color, ligand_color, label,
                          label_offset_mode="camera_up"):
    """Shared rendering logic for both figures - camera locks onto the actual
    M790<->ligand contact zone instead of the whole pocket, and every element
    gets a visually distinct color.

    label_offset_mode controls how the distance label is positioned:
      - "camera_up": offset computed from the camera's actual up direction
        after the final view is set. Adapts per-structure.
      - "fixed": offset is a static (0, 2.5, 0) vector in the molecule's own
        coordinate frame, independent of camera orientation.
    """
    print(f"[+] Rendering {label} using {pdb_file.name}...")
    setup_pymol_display()

    cmd.load(str(pdb_file), "complex")
    obj = "complex"

    resn_selector = " or ".join(f"resn {r}" for r in ligand_resn_list)
    cmd.select("ligand", f"{obj} and ({resn_selector})")
    cmd.select("m790", f"{obj} and resi 790")

    if cmd.count_atoms("ligand") == 0:
        print(f"  [!] Warning: no ligand atoms matched {ligand_resn_list} - check residue naming in the PDB.")
    if cmd.count_atoms("m790") == 0:
        print("  [!] Warning: no atoms matched resi 790 - check residue numbering in the PDB.")

    # Measure the actual closest heavy-atom distance between M790 and the
    # ligand - this is the real number behind a "clash" or "fit" claim, not
    # a guess based on a loose 5 A cutoff.
    result = closest_atom_pair("m790", "ligand")
    if result is not None:
        min_dist, atom_m790, atom_ligand = result
        print(f"  [i] Closest heavy-atom contact: M790/{atom_m790.name} <-> "
              f"{atom_ligand.resn}/{atom_ligand.name}  =  {min_dist:.2f} A")
        if min_dist < 3.2:
            print("  [i] This is tight enough to read as a genuine steric clash.")
        elif min_dist < 4.2:
            print("  [i] This is a close contact, but surfaces will show a small gap, not overlap.")
        else:
            print("  [i] This is NOT a close contact - the residue and ligand are not "
                  "sterically interacting at this distance. Consider whether the "
                  "'collision'/'accommodation' framing is supported by this structure.")

        # Build a tight selection around that exact closest pair (with a small
        # radius so we get a bit of surrounding atoms for context, not just
        # two lone points) and center the camera on THAT, not a loose 5 A shell.
        cmd.select("closest_m790", f"m790 and id {atom_m790.id}")
        cmd.select("closest_ligand", f"ligand and id {atom_ligand.id}")
        cmd.select("contact_m790", "m790 within 4.5 of ligand")
        cmd.select("contact_ligand", "ligand within 4.5 of m790")
        contact_atoms = cmd.count_atoms("contact_m790") + cmd.count_atoms("contact_ligand")
    else:
        min_dist = None
        contact_atoms = 0
        print("  [!] Could not measure distance - m790 or ligand selection is empty.")

    # Background context residues, kept visually minimal
    cmd.select("pocket", f"byres ({obj} within 5.0 of (ligand or m790)) and not m790")

    cmd.hide("everything")

    # Cartoon: distinct blue, very transparent so it reads as background only
    cmd.show("cartoon", obj)
    cmd.color(COLOR_CARTOON, obj)
    cmd.set("cartoon_transparency", 0.75, obj)

    # Pocket context residues: distinct yellow, thin sticks, no surface
    cmd.show("sticks", "pocket")
    cmd.color(COLOR_POCKET, "pocket and elem C")
    cmd.set("stick_radius", 0.10, "pocket")

    # M790: solid color + solid-ish surface so its actual volume is visible
    cmd.show("sticks", "m790")
    cmd.color(m790_color, "m790 and elem C")
    cmd.show("surface", "m790")
    cmd.color(m790_color, "m790")
    cmd.set("transparency", 0.2, "m790")

    # Ligand: solid color + surface, so overlap (collision) or gap
    # (accommodation) against the M790 surface is visually obvious
    cmd.show("sticks", "ligand")
    cmd.color(ligand_color, "ligand and elem C")
    cmd.set("stick_radius", 0.24, "ligand")
    cmd.show("surface", "ligand")
    cmd.color(ligand_color, "ligand")
    cmd.set("transparency", 0.45, "ligand")

    # Draw the actual measured distance as a dashed line with an Angstrom
    # label directly on the figure. For a modest-but-real difference like
    # 3.67 A vs 4.31 A, this is more honest and more convincing than
    # relying on surface overlap alone - the reader sees the real number.
    # NOTE: this must run AFTER hide("everything") / the show() calls above,
    # or hide("everything") wipes the dash/label out before it ever renders.
    if min_dist is not None:
        cmd.distance("contact_dist", "closest_m790", "closest_ligand")
        cmd.hide("labels", "contact_dist")  # hide default atom-pair label
        cmd.color("black", "contact_dist")
        cmd.set("dash_width", 4, "contact_dist")
        cmd.set("dash_gap", 0.3, "contact_dist")
        cmd.set("dash_radius", 0.04, "contact_dist")

    # Camera: center precisely on the single closest atom pair (the actual
    # point of contact/near-contact) rather than a loose shell, so the
    # figure isn't just "in the neighborhood" of the interaction.
    if min_dist is not None:
        cmd.orient("closest_m790 or closest_ligand")
        cmd.zoom("closest_m790 or closest_ligand", buffer=6.0)
    elif contact_atoms > 0:
        cmd.orient("contact_m790 or contact_ligand")
        cmd.zoom("contact_m790 or contact_ligand", buffer=3.0)
    else:
        print("  [!] No contact atoms found - falling back to a wider view. "
              "Check chain/residue naming.")
        cmd.orient("m790 or ligand")
        cmd.zoom("m790 or ligand", buffer=4.0)

    cmd.turn("y", 15)

    # Place the distance label. Two modes are supported (see docstring):
    # "camera_up" adapts to each structure's final camera orientation;
    # "fixed" uses a static molecule-space offset regardless of camera.
    if min_dist is not None:
        ligand_coord = np.array(cmd.get_atom_coords("closest_ligand"))

        if label_offset_mode == "fixed":
            anchor_pos = ligand_coord + np.array([0.0, 2.5, 0.0])
        else:  # "camera_up"
            view = cmd.get_view()
            rot = view[0:9]  # 3x3 model->camera rotation, row-major
            # Model-space direction that maps to camera-space "up" (0,1,0) is
            # the second COLUMN of this row-major matrix: (rot[1], rot[4], rot[7])
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
        cmd.set("label_font_id", 7)  # bold sans-serif, readable at print size
        cmd.label("dist_label_anchor", f'"{min_dist:.2f} A"')
        cmd.show("labels", "dist_label_anchor")

    cmd.ray(2400, 1800)
    cmd.png(str(output_png), dpi=300)
    print(f"  [Saved] {label} -> {output_png}")


def render_figure_2a(erlotinib_pdb, output_png):
    render_pocket_figure(
        erlotinib_pdb, output_png,
        ligand_resn_list=["UNK", "LIG", "ERL"],
        m790_color=COLOR_M790_A,
        ligand_color=COLOR_LIGAND_A,
        label="Figure 2A (Erlotinib collision)",
        label_offset_mode="fixed",  # reverted per feedback: this placement was clean for this structure
    )


def render_figure_2b(osimertinib_pdb, output_png):
    render_pocket_figure(
        osimertinib_pdb, output_png,
        ligand_resn_list=["UNK", "LIG", "OSI"],
        m790_color=COLOR_M790_B,
        ligand_color=COLOR_LIGAND_B,
        label="Figure 2B (Osimertinib accommodation)",
        label_offset_mode="camera_up",  # kept as-is - working well
    )


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    try:
        erlotinib_pdb = find_pdb_file("erlotinib")
        osimertinib_pdb = find_pdb_file("osimertinib")
    except FileNotFoundError as err:
        print(f"[-] File Resolution Error: {err}")
        sys.exit(1)

    fig2a_path = OUTPUT_DIR / "Figure2A_Erlotinib_M790_Collision.png"
    fig2b_path = OUTPUT_DIR / "Figure2B_Osimertinib_M790_Accommodation.png"

    render_figure_2a(erlotinib_pdb, fig2a_path)
    render_figure_2b(osimertinib_pdb, fig2b_path)

    print("\n[+] Both PyMOL figures successfully rendered and exported to results/figures/")


if __name__ == "__main__":
    main()