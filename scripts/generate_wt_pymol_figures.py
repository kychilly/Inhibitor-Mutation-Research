#!/usr/bin/env python3
"""
Automates PyMOL rendering for manuscript wildtype figures with high-visibility
interaction lines and post-rendered distance labels placed directly in front of all 3D geometry.
"""

import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

try:
    import pymol
    from pymol import cmd
except ImportError:
    print("[-] Error: PyMOL module not found. Please run inside an environment with PyMOL installed.")
    sys.exit(1)

OUTPUT_DIR = Path("results/figures")
PLIP_DIR = Path("results/plip_results/wildtype")

COLOR_CARTOON = "skyblue"
COLOR_POCKET = "yellow"
COLOR_T790_A = "deepteal"
COLOR_T790_B = "forest"
COLOR_LIGAND_A = "cyan"
COLOR_LIGAND_B = "orange"

# High-contrast color dedicated specifically to interaction measurements
COLOR_MEASURE = (255, 20, 147)  # HotPink (RGB)


def find_pdb_file(drug_name):
    drug_dir = PLIP_DIR / drug_name
    if drug_dir.exists():
        candidates = sorted(list(drug_dir.glob(f"plipfixed.*{drug_name}_complex_*.pdb")))
        if not candidates:
            candidates = sorted(list(drug_dir.glob(f"*{drug_name}_complex*.pdb")))
    else:
        candidates = sorted(list(PLIP_DIR.glob(f"plipfixed.*{drug_name}_complex_*.pdb")))
        if not candidates:
            candidates = sorted(list(PLIP_DIR.glob(f"*{drug_name}_complex*.pdb")))

    if not candidates:
        raise FileNotFoundError(f"Could not locate Wildtype PDB file for '{drug_name}' in {PLIP_DIR}")

    chosen = candidates[0]
    print(f"  [+] Found PDB for {drug_name}: {chosen.relative_to(PLIP_DIR)}")
    return chosen


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


def render_wt_pocket_figure(pdb_file, output_png, ligand_resn_list, t790_color, ligand_color, label):
    """Shared rendering logic with high-visibility distance dash and 2D post-rendered text."""
    print(f"[+] Rendering {label} using {pdb_file.name}...")
    setup_pymol_display()

    cmd.load(str(pdb_file), "complex")
    obj = "complex"

    resn_selector = " or ".join(f"resn {r}" for r in ligand_resn_list)
    cmd.select("ligand", f"{obj} and ({resn_selector})")
    cmd.select("t790", f"{obj} and resi 790")

    result = closest_atom_pair("t790", "ligand")
    if result is not None:
        min_dist, atom_t790, atom_ligand = result
        cmd.select("closest_t790", f"t790 and id {atom_t790.id}")
        cmd.select("closest_ligand", f"ligand and id {atom_ligand.id}")
        cmd.select("contact_t790", "t790 within 4.5 of ligand")
        cmd.select("contact_ligand", "ligand within 4.5 of t790")
    else:
        min_dist = None

    cmd.select("pocket", f"byres ({obj} within 5.0 of (ligand or t790)) and not t790")
    cmd.hide("everything")

    # Cartoon backbone
    cmd.show("cartoon", obj)
    cmd.color(COLOR_CARTOON, obj)
    cmd.set("cartoon_transparency", 0.75, obj)

    # Context pocket sticks
    cmd.show("sticks", "pocket")
    cmd.color(COLOR_POCKET, "pocket and elem C")
    cmd.set("stick_radius", 0.10, "pocket")

    # T790 residue
    cmd.show("sticks", "t790")
    cmd.color(t790_color, "t790 and elem C")
    cmd.show("surface", "t790")
    cmd.color(t790_color, "t790")
    cmd.set("transparency", 0.2, "t790")

    # Ligand
    cmd.show("sticks", "ligand")
    cmd.color(ligand_color, "ligand and elem C")
    cmd.set("stick_radius", 0.24, "ligand")
    cmd.show("surface", "ligand")
    cmd.color(ligand_color, "ligand")
    cmd.set("transparency", 0.45, "ligand")

    # High-contrast distance dash line
    if min_dist is not None:
        cmd.distance("contact_dist", "closest_t790", "closest_ligand")
        cmd.hide("labels", "contact_dist")  # Hide built-in distance text
        cmd.color("hotpink", "contact_dist")
        cmd.set("dash_width", 5, "contact_dist")
        cmd.set("dash_gap", 0.25, "contact_dist")
        cmd.set("dash_radius", 0.05, "contact_dist")

    # Camera Framing
    if min_dist is not None:
        cmd.orient("closest_t790 or closest_ligand")
        cmd.zoom("closest_t790 or closest_ligand", buffer=6.0)
    else:
        cmd.orient("t790 or ligand")
        cmd.zoom("t790 or ligand", buffer=4.0)

    cmd.turn("y", 15)

    width, height = 2400, 1800
    cmd.ray(width, height)
    cmd.png(str(output_png), dpi=300)

    # Compute exact midpoint in pixel space and overlay distance string over the final PNG
    if min_dist is not None:
        ligand_coord = np.array(cmd.get_atom_coords("closest_ligand"))
        t790_coord = np.array(cmd.get_atom_coords("closest_t790"))
        midpoint = (ligand_coord + t790_coord) / 2.0

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

    wt_erl_path = OUTPUT_DIR / "Figure_WT_Erlotinib_Clearance.png"
    wt_osi_path = OUTPUT_DIR / "Figure_WT_Osimertinib_Clearance.png"

    render_wt_pocket_figure(
        erlotinib_pdb, wt_erl_path,
        ligand_resn_list=["UNK", "LIG", "ERL"],
        t790_color=COLOR_T790_A,
        ligand_color=COLOR_LIGAND_A,
        label="Figure WT-Erlotinib",
    )

    render_wt_pocket_figure(
        osimertinib_pdb, wt_osi_path,
        ligand_resn_list=["UNK", "LIG", "OSI"],
        t790_color=COLOR_T790_B,
        ligand_color=COLOR_LIGAND_B,
        label="Figure WT-Osimertinib",
    )

    print("\n[+] Both wildtype figures rendered with 2D post-overlay measurement labels.")


if __name__ == "__main__":
    main()