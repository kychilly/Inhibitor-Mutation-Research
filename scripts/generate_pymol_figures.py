#!/usr/bin/env python3
# This is an abomination pls help
"""
Automates PyMOL rendering for manuscript figures:
  - Figure 2A: T790M mutant complexed with Erlotinib showing steric collision
               between the bulky Methionine 790 side-chain and Erlotinib's aniline ring.
  - Figure 2B: T790M mutant complexed with Osimertinib showing structural
               accommodation of the M790 side chain by Osimertinib's flexible scaffold.
"""

import sys
from pathlib import Path

# Ensure PyMOL Python module is available
try:
    import pymol
    from pymol import cmd
except ImportError:
    print("[-] Error: PyMOL module not found. Please run this inside an environment with PyMOL installed.")
    sys.exit(1)

OUTPUT_DIR = Path("results/figures")
PLIP_DIR = Path("results/plip_results/t790m")


def find_pdb_file(drug_name):
    """Locates the intermediate plipfixed PDB file for T790M and a given drug."""
    drug_dir = PLIP_DIR / drug_name
    if not drug_dir.exists():
        candidates = list(PLIP_DIR.glob(f"plipfixed.t790m_{drug_name}_complex_*.pdb"))
    else:
        candidates = list(drug_dir.glob(f"plipfixed.t790m_{drug_name}_complex_*.pdb"))

    if not candidates:
        raise FileNotFoundError(f"Could not locate PDB file for T790M + {drug_name} in {PLIP_DIR}")

    return candidates[0]


def setup_pymol_display():
    """Configures global PyMOL aesthetic parameters for camera-ready rendering."""
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 0)          # shadows tend to muddy tight close-ups
    cmd.set("antialias", 2)
    cmd.set("stick_radius", 0.20)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.set("surface_quality", 2)
    cmd.set("transparency_mode", 2)    # per-object transparency blends more predictably
    cmd.set("depth_cue", 0)            # avoid fog washing out a tight close-up
    cmd.set_color("darkorange", [1.0, 0.55, 0.0])  # not a built-in PyMOL color name


def render_figure_2a(erlotinib_pdb, output_png):
    """Figure 2A: M790 side chain visibly clashing into Erlotinib's aniline ring."""
    print(f"[+] Rendering Figure 2A using {erlotinib_pdb.name}...")
    setup_pymol_display()

    cmd.load(str(erlotinib_pdb), "t790m_erlotinib")
    obj = "t790m_erlotinib"

    cmd.select("ligand", f"{obj} and (resn UNK or resn LIG or resn ERL)")
    cmd.select("m790", f"{obj} and resi 790")
    # Just the handful of residues immediately ringing the contact zone -
    # enough for context, not enough to bury the clash visually.
    cmd.select("pocket", f"byres ({obj} within 5.0 of (ligand or m790)) and not m790")

    cmd.hide("everything")

    # Faint cartoon for orientation only - kept very light so it doesn't compete
    cmd.show("cartoon", obj)
    cmd.color("gray90", obj)
    cmd.set("cartoon_transparency", 0.7, obj)

    # Nearby pocket residues: thin, muted, background context
    cmd.show("sticks", "pocket")
    cmd.color("gray70", "pocket and elem C")
    cmd.set("stick_radius", 0.12, "pocket")

    # M790: the star of the clash - solid sticks + solid surface
    cmd.show("sticks", "m790")
    cmd.color("firebrick", "m790 and elem C")
    cmd.show("surface", "m790")
    cmd.color("firebrick", "m790")
    cmd.set("transparency", 0.25, "m790")

    # Ligand: sticks + its own surface, so the two surfaces visibly overlap/press together
    cmd.show("sticks", "ligand")
    cmd.color("teal", "ligand and elem C")
    cmd.set("stick_radius", 0.24, "ligand")
    cmd.show("surface", "ligand")
    cmd.color("teal", "ligand")
    cmd.set("transparency", 0.45, "ligand")

    # Tight framing directly on the clash interface, not the whole pocket
    cmd.orient("m790 or ligand")
    cmd.zoom("m790 or ligand", buffer=1.5)
    cmd.turn("y", 15)  # small offset angle so both surfaces read as 3D, not flat overlap

    cmd.ray(2400, 1800)
    cmd.png(str(output_png), dpi=300)
    print(f"  [Saved] Figure 2A -> {output_png}")


def render_figure_2b(osimertinib_pdb, output_png):
    """Figure 2B: M790 comfortably accommodated by Osimertinib's flexible scaffold."""
    print(f"[+] Rendering Figure 2B using {osimertinib_pdb.name}...")
    setup_pymol_display()

    cmd.load(str(osimertinib_pdb), "t790m_osimertinib")
    obj = "t790m_osimertinib"

    cmd.select("ligand", f"{obj} and (resn UNK or resn LIG or resn OSI)")
    cmd.select("m790", f"{obj} and resi 790")
    cmd.select("pocket", f"byres ({obj} within 5.0 of (ligand or m790)) and not m790")

    cmd.hide("everything")

    cmd.show("cartoon", obj)
    cmd.color("gray90", obj)
    cmd.set("cartoon_transparency", 0.7, obj)

    cmd.show("sticks", "pocket")
    cmd.color("gray70", "pocket and elem C")
    cmd.set("stick_radius", 0.12, "pocket")

    # M790: same solid treatment as 2A so the two figures are visually comparable
    cmd.show("sticks", "m790")
    cmd.color("forest", "m790 and elem C")
    cmd.show("surface", "m790")
    cmd.color("forest", "m790")
    cmd.set("transparency", 0.25, "m790")

    # Ligand surface too - here the point is a visible GAP/fit around the scaffold,
    # not overlap, so keep it a bit more transparent to let the gap read clearly
    cmd.show("sticks", "ligand")
    cmd.color("darkorange", "ligand and elem C")
    cmd.set("stick_radius", 0.24, "ligand")
    cmd.show("surface", "ligand")
    cmd.color("darkorange", "ligand")
    cmd.set("transparency", 0.55, "ligand")

    cmd.orient("m790 or ligand")
    cmd.zoom("m790 or ligand", buffer=1.5)
    cmd.turn("y", 15)  # same offset angle as 2A for a matched camera pair

    cmd.ray(2400, 1800)
    cmd.png(str(output_png), dpi=300)
    print(f"  [Saved] Figure 2B -> {output_png}")


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