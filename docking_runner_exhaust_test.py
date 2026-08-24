#!/usr/bin/env python3
"""
docking_runner_exhaust_test.py

Benchmarks AutoDock Vina search space convergence by running Osimertinib
against Wildtype EGFR across varying exhaustiveness levels (2, 4, 8, 16, 32).
"""

import io
import json
import subprocess
import yaml
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO, Select
from openbabel import openbabel as ob
import pandas as pd

CONFIG_PATH = Path("config/docking_config.yaml")
VINA_BIN = Path("bin/vina.exe")

PROCESSED_DIR = Path("data/processed")
RECEPTOR_DIR = PROCESSED_DIR / "receptors"
LIGAND_DIR = PROCESSED_DIR / "ligands"
RESULTS_DIR = PROCESSED_DIR / "exhaustiveness_test_results"


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

    rec_key = "wildtype"
    lig_name = "osimertinib"
    exhaustiveness_levels = [2, 4, 8, 16, 32]
    fixed_seed = 1

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Phase 2B: Search Space Exhaustiveness Convergence Benchmark ===")

    if rec_key not in receptors_cfg:
        print(f"[-] Receptor key '{rec_key}' not found in configuration.")
        return

    rec_info = receptors_cfg[rec_key]
    rec_path = Path(rec_info["prepared_path"])

    if not rec_path.exists():
        raw_cif = Path(rec_info["raw_path"])
        if raw_cif.exists():
            print(f"[!] Prepared receptor {rec_path} missing. Converting from {raw_cif}...")
            rec_path.parent.mkdir(parents=True, exist_ok=True)
            convert_cif_to_pdbqt_in_memory(str(raw_cif), str(rec_path))
        else:
            print(f"[-] Receptor raw file missing: {raw_cif}")
            return

    lig_path = LIGAND_DIR / f"{lig_name}.pdbqt"
    if not lig_path.exists():
        print(f"[-] Ligand PDBQT missing: {lig_path}")
        return

    center_x = rec_info["center_x"]
    center_y = rec_info["center_y"]
    center_z = rec_info["center_z"]

    print(f"\n[+] Target: {rec_key.upper()} | Ligand: {lig_name.upper()}")
    print(f"[+] Grid Center: ({center_x}, {center_y}, {center_z})")
    print(f"[+] Fixed Seed: {fixed_seed} | Testing Levels: {exhaustiveness_levels}\n")

    results_summary = []

    for exhaust in exhaustiveness_levels:
        out_pdbqt = RESULTS_DIR / f"{rec_key}_{lig_name}_exhaust{exhaust}.pdbqt"
        out_log = RESULTS_DIR / f"{rec_key}_{lig_name}_exhaust{exhaust}.log"

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
            "--exhaustiveness", str(exhaust),
            "--num_modes", str(dock_cfg["num_modes"]),
            "--energy_range", str(dock_cfg["energy_range"]),
            "--seed", str(fixed_seed),
            "--out", str(out_pdbqt)
        ]

        print(f"Running Exhaustiveness = {exhaust:2d} ... ", end="", flush=True)

        with open(out_log, "w") as log_file:
            process = subprocess.run(cmd, stdout=log_file, stderr=subprocess.PIPE, text=True)

        if process.returncode != 0:
            print(f"[FAILED] {process.stderr.strip()}")
            results_summary.append((exhaust, "ERROR"))
        else:
            top_score = parse_top_affinity(out_log)
            print(f"[COMPLETED] Top ΔG: {top_score:.2f} kcal/mol")
            results_summary.append((exhaust, top_score))

    # Save outputs after all levels complete
    csv_out = RESULTS_DIR / "exhaustiveness_convergence.csv"
    json_out = RESULTS_DIR / "exhaustiveness_convergence.json"

    df = pd.DataFrame(results_summary, columns=["exhaustiveness", "top_delta_g"])
    df.to_csv(csv_out, index=False)

    with open(json_out, "w", encoding="utf-8") as f:
        json.dump(dict(results_summary), f, indent=2)

    print(f"\n[+] Saved CSV summary  -> {csv_out}")
    print(f"[+] Saved JSON summary -> {json_out}")

    # Print Summary Convergence Table
    print("\n" + "=" * 45)
    print(" EXHAUSTIVENESS CONVERGENCE SUMMARY TABLE ")
    print("=" * 45)
    print(f"{'Exhaustiveness':<18} | {'Top ΔG (kcal/mol)':<20}")
    print("-" * 45)
    for exhaust, score in results_summary:
        score_str = f"{score:.2f}" if isinstance(score, float) else str(score)
        print(f"{exhaust:<18} | {score_str:<20}")
    print("=" * 45 + "\n")


if __name__ == "__main__":
    main()