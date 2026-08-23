import os
import numpy as np
import yaml
from pathlib import Path
from Bio.PDB import MMCIFParser


def extract_native_ligand_center(cif_path, ligand_resname="AQ4"):
    """
    Parses the 1M17 mmCIF file, locates the native ligand residue (AQ4),
    and calculates its 3D geometric center (centroid).
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("1M17", str(cif_path))

    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Match residue name (1M17 native ligand is AQ4/Erlotinib)
                if residue.get_resname().strip() == ligand_resname:
                    for atom in residue:
                        coords.append(atom.get_coord())

    if not coords:
        raise ValueError(f"Ligand '{ligand_resname}' not found in {cif_path}.")

    coords = np.array(coords)
    centroid = np.mean(coords, axis=0)
    return [round(float(c), 3) for c in centroid]


def main():
    raw_cif = Path("data/raw/receptors/1M17.cif")
    config_path = Path("config/docking_config.yaml")

    if not raw_cif.exists():
        print(f"[-] Error: Raw 1M17 CIF not found at {raw_cif}")
        return

    print(f"[+] Extracting native ligand geometric center from {raw_cif}...")
    centroid = extract_native_ligand_center(raw_cif, ligand_resname="AQ4")

    print(f"[+] Computed Pocket Center (x, y, z): {centroid}")

    # Load existing config if available
    if config_path.exists():
        with open(config_path, "r") as f:
            config = yaml.safe_load(f) or {}
    else:
        config = {}

    # Lock in standard grid box parameters
    config["grid"] = {
        "center_x": centroid[0],
        "center_y": centroid[1],
        "center_z": centroid[2],
        "size_x": 22.0,
        "size_y": 22.0,
        "size_z": 22.0
    }

    # Explicitly set accurate receptor paths to data/raw/receptors/
    config["receptors"] = {
        "wildtype": "data/raw/receptors/1M17.cif",
        "t790m": "data/raw/receptors/2JIT.cif"
    }

    # Write back to config file
    config_path.parent.mkdir(parents=True, exist_ok=True)
    with open(config_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    print(f"[+] Successfully locked grid box window and receptor paths in {config_path}")


if __name__ == "__main__":
    main()