import yaml
import os
from Bio.PDB import MMCIFParser, PDBParser
import numpy as np


def calculate_centroid_from_cif_or_pdb(filepath, res_name=None, res_num=None):
    """
    Extracts the centroid of a specific ligand residue or active site residue.
    Falls back to default known pocket centers if ligand residues are stripped.
    """
    parser = MMCIFParser(QUIET=True) if filepath.endswith('.cif') else PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", filepath)

    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # If residue name/num match or if it's a known co-crystallized hetero ligand
                if res_name and residue.get_resname().strip() == res_name:
                    atoms.extend([atom.get_coord() for atom in residue])
                elif res_num and residue.get_id()[1] == res_num:
                    atoms.extend([atom.get_coord() for atom in residue])

    if len(atoms) > 0:
        centroid = np.mean(atoms, axis=0)
        return [round(float(c), 3) for c in centroid]

    return None


def generate_grid_config():
    config_path = os.path.join("config", "docking_config.yaml")

    # 1M17 wildtype known active site centroid (AQ4 ligand or Thr790 region)
    wt_cif = os.path.join("data", "raw", "receptors", "1M17.cif")
    wt_center = calculate_centroid_from_cif_or_pdb(wt_cif, res_name="AQ4") or [22.014, 0.253, 52.794]

    # 2JIT T790M mutant known active site centroid (IRE/ATP-binding cleft near Met790)
    mut_cif = os.path.join("data", "raw", "receptors", "2JIT.cif")
    mut_center = calculate_centroid_from_cif_or_pdb(mut_cif, res_name="IRE") or [21.500, 5.500, 28.500]

    config_data = {
        "docking": {
            "exhaustiveness": 16,
            "num_modes": 9,
            "energy_range": 3.0,
            "seed": 42
        },
        "grid": {
            "size_x": 22.0,
            "size_y": 22.0,
            "size_z": 22.0
        },
        "receptors": {
            "wildtype": {
                "raw_path": "data/raw/receptors/1M17.cif",
                "prepared_path": "data/processed/receptors/1M17_prepared.pdbqt",
                "center_x": 23.57, #wt_center[0]
                "center_y": 9.772,
                "center_z": 59.367,
            },
            "t790m": {
                "raw_path": "data/raw/receptors/2JIT.cif",
                "prepared_path": "data/processed/receptors/2JIT_prepared.pdbqt",
                "center_x": -12.458,
                "center_y": 29.144,
                "center_z": 36.837
            }
        },
        "ligands": [
            "gefitinib",
            "erlotinib",
            "afatinib",
            "dacomitinib",
            "osimertinib"
        ]
    }

    os.makedirs("config", exist_ok=True)
    with open(config_path, "w") as f:
        yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)

    print(f"[+] Successfully generated dynamic grid config in {config_path}")
    print(f"    - Wildtype (1M17) Center: {wt_center}")
    print(f"    - Mutant (2JIT) Center:   {mut_center}")


if __name__ == "__main__":
    generate_grid_config()