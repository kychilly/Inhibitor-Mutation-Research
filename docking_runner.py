import os
import csv
import subprocess
import yaml
import numpy as np

# Path configurations
CONFIG_PATH = "config/docking_config.yaml"
VINA_BIN = os.path.join("bin", "vina.exe")

# Directories
PROCESSED_DIR = "data/processed"
PROCESSED_RECEPTORS_DIR = os.path.join(PROCESSED_DIR, "receptors")
PROCESSED_LIGANDS_DIR = os.path.join(PROCESSED_DIR, "ligands")

RESULTS_DIR = os.path.join(PROCESSED_DIR, "docking_results")
RAW_DATA_CSV = os.path.join(RESULTS_DIR, "raw_docking_scores.csv")
NUM_SEEDS = 15


def load_config(config_path):
    """Loads parameters from docking_config.yaml."""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def verify_grid_alignment(receptor_path, grid_center, max_dist=6.0):
    """
    SAFEGUARD: Parses PDBQT receptor file for Gatekeeper residue 790
    and calculates distance to grid center.
    """
    if not os.path.exists(receptor_path):
        print(f"  [Safeguard Notice] Receptor '{receptor_path}' not found. Skipping check.")
        return

    gatekeeper_coords = []
    with open(receptor_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                # Flexibly scan for residue 790 across PDBQT column variations
                if len(parts) >= 8:
                    # Check if residue number (typically token 5 or 6) is 790
                    if "790" in parts[4:7]:
                        try:
                            # Extract X, Y, Z coordinates from standard PDB coordinate columns
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            gatekeeper_coords.append([x, y, z])
                        except ValueError:
                            continue

    if not gatekeeper_coords:
        print(f"  [Safeguard Notice] Residue 790 not found in '{receptor_path}'. Proceeding.")
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


def parse_top_affinity(log_file_path):
    """Extracts top pose score (kcal/mol) from AutoDock Vina log output."""
    with open(log_file_path, "r") as log:
        for line in log:
            if line.strip().startswith("1 "):  # Top pose entry starts with mode '1'
                parts = line.split()
                return float(parts[1])
    raise RuntimeError(f"Could not parse top pose affinity score from {log_file_path}")


def resolve_ligand_path(lig_name):
    """Dynamically finds the ligand .pdbqt file across potential directory layouts."""
    possible_paths = [
        os.path.join(PROCESSED_LIGANDS_DIR, f"{lig_name}_minimized.pdbqt"),
        os.path.join(PROCESSED_LIGANDS_DIR, f"{lig_name}.pdbqt"),
        os.path.join(PROCESSED_DIR, f"{lig_name}_minimized.pdbqt"),
        os.path.join(PROCESSED_DIR, f"{lig_name}.pdbqt")
    ]
    for path in possible_paths:
        if os.path.exists(path):
            return path
    raise FileNotFoundError(
        f"Missing processed ligand file for '{lig_name}'. "
        f"Checked: {possible_paths}"
    )


def main():
    if not os.path.exists(VINA_BIN):
        raise FileNotFoundError(
            f"Vina binary not found at '{VINA_BIN}'. Please place vina.exe in the 'bin' directory."
        )

    os.makedirs(RESULTS_DIR, exist_ok=True)
    cfg = load_config(CONFIG_PATH)

    grid = cfg["grid"]
    docking_opts = cfg["docking"]
    ligand_names = cfg["ligands"]

    # Map receptor keys to PROCESSED PDBQT paths
    receptor_paths = {
        "wildtype": os.path.join(PROCESSED_RECEPTORS_DIR, "1M17_prepared.pdbqt"),
        "t790m": os.path.join(PROCESSED_RECEPTORS_DIR, "2JIT_prepared.pdbqt")
    }

    print("=== Phase 1: Running Grid Safeguard Checks ===")
    for rec_type, rec_path in receptor_paths.items():
        verify_grid_alignment(rec_path, grid)

    print("\n=== Phase 2: Executing Docking Matrix ===")
    raw_results = []

    for rec_key, rec_path in receptor_paths.items():
        if not os.path.exists(rec_path):
            raise FileNotFoundError(f"Missing processed receptor file: '{rec_path}'")

        for lig_name in ligand_names:
            lig_path = resolve_ligand_path(lig_name)

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
                    "--exhaustiveness", str(docking_opts.get("exhaustiveness", 32)),
                    "--num_modes", str(docking_opts.get("num_modes", 9)),
                    "--energy_range", str(docking_opts.get("energy_range", 3.0)),
                    "--seed", str(run_seed),
                    "--out", out_pdbqt,
                    "--log", log_file
                ]

                subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
                top_score = parse_top_affinity(log_file)

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