import os
import csv
import subprocess
import yaml
import io
import numpy as np
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO, Select
from openbabel import openbabel as ob

CONFIG_PATH = Path("config/docking_config.yaml")
VINA_BIN = Path("bin/vina.exe")

PROCESSED_DIR = Path("data/processed")
RECEPTOR_DIR = PROCESSED_DIR / "receptors"
LIGAND_DIR = PROCESSED_DIR / "ligands"
RESULTS_DIR = PROCESSED_DIR / "docking_results"

class CleanProteinSelect(Select):
    """Retains standard protein residues while stripping water and heteroatoms."""
    def accept_residue(self, residue):
        return residue.get_resname() not in ["HOH", "WAT"] and residue.id[0] == " "

def sanitize_rigid_receptor_pdbqt(pdbqt_string):
    forbidden_tags = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")
    cleaned_lines = []
    for line in pdbqt_string.splitlines():
        if not line.strip().startswith(forbidden_tags):
            cleaned_lines.append(line)
    return "\n".join(cleaned_lines) + "\n"

def convert_cif_to_pdbqt_in_memory(cif_path, pdbqt_path):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", cif_path)

    pdb_io = PDBIO()
    pdb_io.set_structure(structure)

    string_stream = io.StringIO()
    pdb_io.save(string_stream, select=CleanProteinSelect())
    pdb_string = string_stream.getvalue()

    ob_conversion = ob.OBConversion()
    ob_conversion.SetInAndOutFormats("pdb", "pdbqt")

    mol = ob.OBMol()
    if not ob_conversion.ReadString(mol, pdb_string):
        raise RuntimeError(f"OpenBabel failed to read structure from {cif_path}")

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

def parse_top_affinity(log_file_path):
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
    cfg = load_config(CONFIG_PATH)
    grid_cfg = cfg["grid"]
    dock_cfg = cfg["docking"]
    receptors_cfg = cfg["receptors"]
    ligands = cfg["ligands"]
    num_seeds = dock_cfg.get("seed", 15) if isinstance(dock_cfg.get("seed"), int) else 15

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Phase 2: Executing Docking Matrix ===")

    # Run for t790m (or iterate through receptors_cfg for both)
    target_receptors = ["t790m"]  # Change to ["wildtype", "t790m"] after this single experiment

    for rec_key in target_receptors:
        if rec_key not in receptors_cfg:
            print(f"[-] Receptor key '{rec_key}' not found in configuration.")
            continue

        rec_info = receptors_cfg[rec_key]
        rec_path = Path(rec_info["prepared_path"])

        # Automatically build PDBQT from CIF if missing
        if not rec_path.exists():
            raw_cif = Path(rec_info["raw_path"])
            if raw_cif.exists():
                print(f"[!] Prepared receptor {rec_path} missing. Converting from {raw_cif}...")
                rec_path.parent.mkdir(parents=True, exist_ok=True)
                convert_cif_to_pdbqt_in_memory(str(raw_cif), str(rec_path))
            else:
                print(f"[-] Receptor raw file missing: {raw_cif}")
                continue

        # Extract receptor-specific coordinates from YAML
        center_x = rec_info["center_x"]
        center_y = rec_info["center_y"]
        center_z = rec_info["center_z"]

        print(f"\n[+] Loaded Target: {rec_key.upper()} | Center: ({center_x}, {center_y}, {center_z})")

        for lig_name in ligands:
            lig_path = LIGAND_DIR / f"{lig_name}.pdbqt"
            if not lig_path.exists():
                print(f"[-] Ligand PDBQT missing: {lig_path}")
                continue

            print(f"\nDocking Pair: Receptor='{rec_key}' | Ligand='{lig_name}' across 15 seeds...")

            for seed in range(1, 16):
                out_pdbqt = RESULTS_DIR / f"{rec_key}_{lig_name}_seed{seed}.pdbqt"
                out_log = RESULTS_DIR / f"{rec_key}_{lig_name}_seed{seed}.log"

                cmd = [
                    str(VINA_BIN),
                    "--receptor", str(rec_path),
                    "--ligand", str(lig_path),
                    "--center_x", str(center_x),
                    "--center_y", str(center_y),
                    "--center_z", str(center_z),
                    "--size_x", str(grid_cfg["size_x"]),
                    "--size_y", str(grid_cfg["size_y"]),
                    "--size_z", str(grid_cfg["size_z"]),
                    "--exhaustiveness", str(dock_cfg["exhaustiveness"]),
                    "--num_modes", str(dock_cfg["num_modes"]),
                    "--energy_range", str(dock_cfg["energy_range"]),
                    "--seed", str(seed),
                    "--out", str(out_pdbqt)
                ]

                with open(out_log, "w") as log_file:
                    process = subprocess.run(cmd, stdout=log_file, stderr=subprocess.PIPE, text=True)

                if process.returncode != 0:
                    print(f"  [Error] Seed {seed}/15 failed: {process.stderr.strip()}")
                else:
                    top_score = parse_top_affinity(out_log)
                    score_str = f"{top_score} kcal/mol" if top_score is not None else "N/A"
                    print(f"  [Completed] Seed {seed}/15 | Top Affinity: {score_str}")

if __name__ == "__main__":
    main()