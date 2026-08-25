#!/usr/bin/env python3
"""
Automates PyMOL rendering for manuscript figures:
  - Figure 2A: T790M mutant complexed with Erlotinib showing steric collision
               between the bulky Methionine 790 side-chain and Erlotinib's aniline ring.
  - Figure 2B: T790M mutant complexed with Osimertinib showing structural
               accommodation of the M790 side chain by Osimertinib's flexible scaffold.
"""

import sys
from pathlib import Path

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


def render_pocket_figure(pdb_file, output_png, ligand_resn_list, m790_color, ligand_color, label):
    """Shared rendering logic for both figures - camera locks onto the actual
    M790<->ligand contact zone instead of the whole pocket, and every element
    gets a visually distinct color."""
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

    # The actual contact zone: only the atoms on each side that are close to
    # the other. This is what the camera should center on, not the whole
    # ligand or the whole residue.
    cmd.select("contact_m790", "m790 within 5.0 of ligand")
    cmd.select("contact_ligand", "ligand within 5.0 of m790")
    contact_atoms = cmd.count_atoms("contact_m790") + cmd.count_atoms("contact_ligand")

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

    # Camera: center on the actual contact zone if we found one, otherwise
    # fall back to centering on m790/ligand as a whole so we never render
    # an empty/blank view.
    if contact_atoms > 0:
        cmd.orient("contact_m790 or contact_ligand")
        cmd.zoom("contact_m790 or contact_ligand", buffer=3.0)
    else:
        print("  [!] No atoms found within 5.0 A between M790 and ligand - "
              "falling back to a wider view. Check chain/residue naming.")
        cmd.orient("m790 or ligand")
        cmd.zoom("m790 or ligand", buffer=4.0)

    cmd.turn("y", 15)

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
    )


def render_figure_2b(osimertinib_pdb, output_png):
    render_pocket_figure(
        osimertinib_pdb, output_png,
        ligand_resn_list=["UNK", "LIG", "OSI"],
        m790_color=COLOR_M790_B,
        ligand_color=COLOR_LIGAND_B,
        label="Figure 2B (Osimertinib accommodation)",
    )


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    try:
        erlotinib_pdb = find_pdb_file("erlotinib")
        osimertinib_pdb = find_pdb_file("osimertinib")
    except FileNotFoundError as err:
        print(f"[-] File Resolution Error: {err}")
        sys.exit(1)

    fig2a_path = OUTPUT_DIR / "Figure2A_Erlotinib_M790_Collision_zoomed.png"
    fig2b_path = OUTPUT_DIR / "Figure2B_Osimertinib_M790_Accommodation_zoomed.png"

    render_figure_2a(erlotinib_pdb, fig2a_path)
    render_figure_2b(osimertinib_pdb, fig2b_path)

    print("\n[+] Both PyMOL figures successfully rendered and exported to results/figures/")


if __name__ == "__main__":
    main()