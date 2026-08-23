import os
import csv
import subprocess
import yaml
import io
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
from openbabel import openbabel as ob

CONFIG_PATH = "config/docking_config.yaml"
VINA_BIN = os.path.join("bin", "vina.exe")

PROCESSED_DIR = "data/processed"
PROCESSED_RECEPTORS_DIR = os.path.join(PROCESSED_DIR, "receptors")
PROCESSED_LIGANDS_DIR = os.path.join(PROCESSED_DIR, "ligands")

RESULTS_DIR = os.path.join(PROCESSED_DIR, "docking_results")
RAW_DATA_CSV = os.path.join(RESULTS_DIR, "raw_docking_scores.csv")
NUM_SEEDS = 15


class CleanProteinSelect(Select):
    """Retains standard protein residues while stripping water and heteroatoms."""

    def accept_residue(self, residue):
        return residue.get_resname() not in ["HOH", "WAT"] and residue.id[0] == " "


def sanitize_rigid_receptor_pdbqt(pdbqt_string):
    """
    SAFEGUARD: AutoDock Vina strictly forbids flexibility keywords (ROOT, BRANCH, TORSDOF)
    in rigid receptor PDBQT files.
    """
    forbidden_tags = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")
    cleaned_lines = []
    for line in pdbqt_string.splitlines():
        if not line.strip().startswith(forbidden_tags):
            cleaned_lines.append(line)
    return "\n".join(cleaned_lines) + "\n"


def convert_cif_to_pdbqt_in_memory(cif_path, pdbqt_path):
    """Converts CIF to rigid PDBQT in-memory using io.StringIO."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", cif_path)

    pdb_io = PDBIO()
    pdb_io.set_structure(structure)

    # Standard io.StringIO provides .tell() and .seek() required by Biopython
    string_stream = io.StringIO()
    pdb_io.save(string_stream, select=CleanProteinSelect())
    pdb_string = string_stream.getvalue()

    ob_conversion = ob.OBConversion()
    ob_conversion.SetInAndOutFormats("pdb", "pdbqt")

    mol = ob.OBMol()
    if not ob_conversion.ReadString(mol, pdb_string):
        raise RuntimeError(f"OpenBabel failed to read structure from {cif_path}")

    # Add polar hydrogens at physiological pH (7.4) & calculate Gasteiger partial charges
    mol.AddHydrogens(False, True, 7.4)

    raw_pdbqt_str = ob_conversion.WriteString(mol)
    if not raw_pdbqt_str:
        raise RuntimeError(f"OpenBabel failed to output PDBQT string for {cif_path}")

    clean_pdbqt_str = sanitize_rigid_receptor_pdbqt(raw_pdbqt_str)

    with open(pdbqt_path, "w") as f:
        f.write(clean_pdbqt_str)

    print(f"  [Auto-Prepared Rigid Receptor] {cif_path} -> {pdbqt_path}")


def load_config(config_path):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def verify_grid_alignment(receptor_path, grid_center, max_dist=15.0):
    """
    SAFEGUARD: Parses PDBQT receptor for Gatekeeper residue 790 (T790/M790)
    and verifies that the grid center is within range of the active site.
    """
    if not os.path.exists(receptor_path):
        print(f"  [Safeguard Notice] Receptor '{receptor_path}' not found. Skipping check.")
        return

    gatekeeper_coords = []
    with open(receptor_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                res_seq = line[22:26].strip()
                if res_seq == "790" or "790" in parts[4:7]:
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        gatekeeper_coords.append([x, y, z])
                    except (ValueError, IndexError):
                        continue

    if not gatekeeper_coords:
        print(f"  [Safeguard Warning] Residue 790 not found in '{receptor_path}'. Proceeding with caution.")
        return

    res_centroid = np.mean(gatekeeper_coords, axis=0)
    grid_point = np.array([grid_center["center_x"], grid_center["center_y"], grid_center["center_z"]])
    dist = np.linalg.norm(res_centroid - grid_point)

    if dist > max_dist:
        raise ValueError(
            f"Safeguard Alert: Grid center is {dist:.2f}Å away from Residue 790 "
            f"in '{receptor_path}' (max allowed: {max_dist}Å)!"
        )
    print(f"  [Safeguard Passed] Grid center is {dist:.2f}Å from Residue 790 in {os.path.basename(receptor_path)}.")

def verify_pdbqt_file(pdbqt_path):
    """SAFEGUARD: Checks that PDBQT contains valid non-empty ATOM records."""
    if not os.path.exists(pdbqt_path) or os.path.getsize(pdbqt_path) == 0:
        raise FileNotFoundError(f"PDBQT file missing or empty: {pdbqt_path}")

    has_atoms = False
    with open(pdbqt_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                has_atoms = True
                break
    if not has_atoms:
        raise ValueError(f"No ATOM or HETATM records found in PDBQT file: {pdbqt_path}")


def parse_top_affinity(log_file_path):
    """Extracts top pose score (kcal/mol) from Vina stdout log file."""
    with open(log_file_path, "r") as log:
        for line in log:
            line_str = line.strip()
            if line_str.startswith("1 ") or (len(line_str.split()) > 1 and line_str.split()[0] == "1"):
                parts = line_str.split()
                try:
                    return float(parts[1])
                except (ValueError, IndexError):
                    continue
    raise RuntimeError(f"Could not parse top pose affinity score from {log_file_path}")


def main():
    if not os.path.exists(VINA_BIN):
        raise FileNotFoundError(
            f"Vina binary not found at '{VINA_BIN}'. Please place vina.exe in the 'bin' directory."
        )

    os.makedirs(PROCESSED_RECEPTORS_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    cfg = load_config(CONFIG_PATH)

    grid = cfg["grid"]
    docking_opts = cfg["docking"]
    ligand_names = cfg["ligands"]

    receptor_targets = {
        "wildtype": {
            "cif": os.path.join(PROCESSED_RECEPTORS_DIR, "wildtype_clean.cif"),
            "pdbqt": os.path.join(PROCESSED_RECEPTORS_DIR, "1M17_prepared.pdbqt")
        },
        "t790m": {
            "cif": os.path.join(PROCESSED_RECEPTORS_DIR, "t790m_clean.cif"),
            "pdbqt": os.path.join(PROCESSED_RECEPTORS_DIR, "2JIT_prepared.pdbqt")
        }
    }

    print("=== Phase 0: Checking and Converting Receptors to Rigid PDBQT ===")
    for rec_key, paths in receptor_targets.items():
        needs_prep = True
        if os.path.exists(paths["pdbqt"]):
            with open(paths["pdbqt"], "r") as f:
                content = f.read()
                if "ROOT" not in content and "BRANCH" not in content:
                    needs_prep = False

        if needs_prep:
            if os.path.exists(paths["cif"]):
                convert_cif_to_pdbqt_in_memory(paths["cif"], paths["pdbqt"])
            else:
                raise FileNotFoundError(f"Missing both PDBQT and CIF files for {rec_key}: '{paths['cif']}'")

    print("\n=== Phase 1: Running Grid Safeguard Checks ===")
    for rec_key, paths in receptor_targets.items():
        verify_pdbqt_file(paths["pdbqt"])
        verify_grid_alignment(paths["pdbqt"], grid)

    print("\n=== Phase 2: Executing Docking Matrix ===")
    raw_results = []

    for rec_key, paths in receptor_targets.items():
        rec_path = paths["pdbqt"]

        for lig_name in ligand_names:
            lig_path = os.path.join(PROCESSED_LIGANDS_DIR, f"{lig_name}.pdbqt")
            verify_pdbqt_file(lig_path)

            print(f"\nDocking Pair: Receptor='{rec_key}' | Ligand='{lig_name}' across {NUM_SEEDS} seeds...")

            for seed_idx in range(1, NUM_SEEDS + 1):
                run_seed = docking_opts.get("seed", 42) + seed_idx
                out_pdbqt = os.path.join(RESULTS_DIR, f"{rec_key}_{lig_name}_seed{seed_idx}.pdbqt")
                log_file = os.path.join(RESULTS_DIR, f"{rec_key}_{lig_name}_seed{seed_idx}.log")

                cmd = [
                    VINA_BIN,
                    "--receptor", rec_path,
                    "--ligand", lig_path,
                    "--center_x", str(grid["center_x"]),
                    "--center_y", str(grid["center_y"]),
                    "--center_z", str(grid["center_z"]),
                    "--size_x", str(grid["size_x"]),
                    "--size_y", str(grid["size_y"]),
                    "--size_z", str(grid["size_z"]),
                    "--exhaustiveness", str(docking_opts.get("exhaustiveness", 16)),
                    "--num_modes", str(docking_opts.get("num_modes", 9)),
                    "--energy_range", str(docking_opts.get("energy_range", 3.0)),
                    "--seed", str(run_seed),
                    "--out", out_pdbqt
                ]

                with open(log_file, "w") as log_f:
                    result = subprocess.run(cmd, stdout=log_f, stderr=subprocess.PIPE, text=True)

                if result.returncode != 0:
                    print(f"\n[!] Vina Error on Seed {seed_idx}:")
                    print(result.stderr)
                    raise RuntimeError(f"Vina exited with code {result.returncode}")

                top_score = parse_top_affinity(log_file)
                print(f"  [Completed] Seed {seed_idx}/{NUM_SEEDS} | Top Affinity: {top_score} kcal/mol")
                raw_results.append({
                    "receptor": rec_key,
                    "ligand": lig_name,
                    "seed": seed_idx,
                    "affinity_kcal_mol": top_score
                })

    print(f"\n=== Phase 3: Writing Raw Scores to Disk ===")
    with open(RAW_DATA_CSV, mode="w", newline="") as csv_file:
        fieldnames = ["receptor", "ligand", "seed", "affinity_kcal_mol"]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(raw_results)

    print(f"Successfully saved {len(raw_results)} docking observations to '{RAW_DATA_CSV}'.")


if __name__ == "__main__":
    main()