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

import numpy as np
from PIL import Image, ImageDraw, ImageFont

try:
    import pymol
    from pymol import cmd
except ImportError:
    print("[-] Error: PyMOL module not found. Please run this inside an environment with PyMOL installed.")
    sys.exit(1)

OUTPUT_DIR = Path("results/figures")
PLIP_DIR = Path("results/plip_results/t790m")

# Distinct, high-contrast palette
COLOR_CARTOON = "skyblue"
COLOR_POCKET = "yellow"
COLOR_M790_A = "red"
COLOR_M790_B = "green"
COLOR_LIGAND_A = "cyan"
COLOR_LIGAND_B = "orange"

# High-contrast color dedicated specifically to interaction measurements
COLOR_MEASURE = (255, 20, 147)  # HotPink (RGB)


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


def get_screen_coords(model_coord, img_width=2400, img_height=1800):
    """Map 3D world coordinate to 2D image pixel space (X, Y)."""
    view = cmd.get_view()
    rot = np.array(view[0:9]).reshape((3, 3))
    pos = np.array(view[9:12])
    origin = np.array(view[12:15])

    # Transform 3D coordinate through camera projection
    pt = model_coord - origin
    pt_cam = np.dot(rot, pt) + pos

    # Scale camera matrix field of view to pixel dimension
    fov = cmd.get("field_of_view")
    fov_rad = np.radians(float(fov) if fov is not None else 25.0)
    scale = (img_height / 2.0) / np.tan(fov_rad / 2.0)

    px = (img_width / 2.0) + (pt_cam[0] * scale / abs(pt_cam[2]))
    py = (img_height / 2.0) - (pt_cam[1] * scale / abs(pt_cam[2]))
    return float(px), float(py)


def add_post_render_label(image_path, text, pixel_pos):
    """Draws text with a clean white outline directly over the top layer of the PNG image."""
    img = Image.open(image_path).convert("RGBA")
    draw = ImageDraw.Draw(img)

    try:
        font = ImageFont.truetype("arial.ttf", 64)
    except IOError:
        font = ImageFont.load_default()

    x, y = pixel_pos

    # Get bounding box to center text exactly on point
    bbox = draw.textbbox((0, 0), text, font=font)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]

    # Offset upward slightly so text hovers over the dash midpoint
    draw_x = x - (text_width / 2.0)
    draw_y = y - text_height - 15

    # Draw white text outline (stroke) for high visibility
    stroke_w = 4
    for dx in range(-stroke_w, stroke_w + 1):
        for dy in range(-stroke_w, stroke_w + 1):
            if dx != 0 or dy != 0:
                draw.text((draw_x + dx, draw_y + dy), text, font=font, fill=(255, 255, 255, 255))

    # Draw primary text layer on top of everything
    draw.text((draw_x, draw_y), text, font=font, fill=COLOR_MEASURE + (255,))

    img.save(image_path)


def render_pocket_figure(pdb_file, output_png, ligand_resn_list, m790_color, ligand_color, label):
    print(f"[+] Rendering {label} using {pdb_file.name}...")
    setup_pymol_display()

    cmd.load(str(pdb_file), "complex")
    obj = "complex"

    resn_selector = " or ".join(f"resn {r}" for r in ligand_resn_list)
    cmd.select("ligand", f"{obj} and ({resn_selector})")
    cmd.select("m790", f"{obj} and resi 790")

    result = closest_atom_pair("m790", "ligand")
    if result is not None:
        min_dist, atom_m790, atom_ligand = result
        print(f"  [i] Closest heavy-atom contact: M790/{atom_m790.name} <-> "
              f"{atom_ligand.resn}/{atom_ligand.name}  =  {min_dist:.2f} A")

        cmd.select("closest_m790", f"m790 and id {atom_m790.id}")
        cmd.select("closest_ligand", f"ligand and id {atom_ligand.id}")

    cmd.select("pocket", f"byres ({obj} within 5.0 of (ligand or m790)) and not m790")

    cmd.hide("everything")

    # Cartoon backbone
    cmd.show("cartoon", obj)
    cmd.color(COLOR_CARTOON, obj)
    cmd.set("cartoon_transparency", 0.75, obj)

    # Context pocket
    cmd.show("sticks", "pocket")
    cmd.color(COLOR_POCKET, "pocket and elem C")
    cmd.set("stick_radius", 0.10, "pocket")

    # M790 residue
    cmd.show("sticks", "m790")
    cmd.color(m790_color, "m790 and elem C")
    cmd.show("surface", "m790")
    cmd.color(m790_color, "m790")
    cmd.set("transparency", 0.2, "m790")

    # Ligand
    cmd.show("sticks", "ligand")
    cmd.color(ligand_color, "ligand and elem C")
    cmd.set("stick_radius", 0.24, "ligand")
    cmd.show("surface", "ligand")
    cmd.color(ligand_color, "ligand")
    cmd.set("transparency", 0.45, "ligand")

    # Distance dash line
    if min_dist is not None:
        cmd.distance("contact_dist", "closest_m790", "closest_ligand")
        cmd.hide("labels", "contact_dist")
        cmd.color("hotpink", "contact_dist")
        cmd.set("dash_width", 5, "contact_dist")
        cmd.set("dash_gap", 0.25, "contact_dist")
        cmd.set("dash_radius", 0.05, "contact_dist")

    # Framing
    if min_dist is not None:
        cmd.orient("closest_m790 or closest_ligand")
        cmd.zoom("closest_m790 or closest_ligand", buffer=6.0)
    else:
        cmd.orient("m790 or ligand")
        cmd.zoom("m790 or ligand", buffer=4.0)

    cmd.turn("y", 15)

    width, height = 2400, 1800
    cmd.ray(width, height)
    cmd.png(str(output_png), dpi=300)

    # Compute exact midpoint in pixel space and overlay distance string over the final PNG
    if min_dist is not None:
        ligand_coord = np.array(cmd.get_atom_coords("closest_ligand"))
        m790_coord = np.array(cmd.get_atom_coords("closest_m790"))
        midpoint = (ligand_coord + m790_coord) / 2.0

        screen_x, screen_y = get_screen_coords(midpoint, width, height)
        add_post_render_label(output_png, f"{min_dist:.2f} Å", (screen_x, screen_y))

    print(f"  [Saved] {label} -> {output_png}")


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

    render_pocket_figure(
        erlotinib_pdb, fig2a_path,
        ligand_resn_list=["UNK", "LIG", "ERL"],
        m790_color=COLOR_M790_A,
        ligand_color=COLOR_LIGAND_A,
        label="Figure 2A (Erlotinib collision)",
    )

    render_pocket_figure(
        osimertinib_pdb, fig2b_path,
        ligand_resn_list=["UNK", "LIG", "OSI"],
        m790_color=COLOR_M790_B,
        ligand_color=COLOR_LIGAND_B,
        label="Figure 2B (Osimertinib accommodation)",
    )

    print("\n[+] Both PyMOL figures successfully rendered and exported to results/figures/")


if __name__ == "__main__":
    main()