import os
import yaml
from pathlib import Path
from Bio.PDB import MMCIFParser, MMCIFIO, Select
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem


class ReceptorCleanSelect(Select):
    """
    Select filter to retain only standard protein amino acids (Chain A),
    removing H2O molecules, crystallization ions, and heteroatoms.
    """

    def accept_chain(self, chain):
        return chain.get_id() == 'A'

    def accept_residue(self, residue):
        # Exclude H_HOH (waters) and heteroatom residues (H_*)
        if residue.id[0] != " ":
            return False
        return True


def clean_cif(input_cif_path, output_cif_path):
    """Parses raw mmCIF structure and strips unwanted water/heteroatoms."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", str(input_cif_path))

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(output_cif_path), ReceptorCleanSelect())
    print(f"[+] Cleaned receptor saved to: {output_cif_path}")


def prepare_ligand_pdbqt(sdf_path, output_pdbqt_path):
    """Converts a raw 3D SDF ligand file to AutoDock PDBQT format using Meeko."""
    supplier = Chem.SDMolSupplier(str(sdf_path))
    mol = next(supplier)
    if mol is None:
        raise ValueError(f"Failed to load ligand from {sdf_path}")

    # Add polar hydrogens for docking scoring
    mol = Chem.AddHs(mol, addCoords=True)

    # Prepare parameterization via Meeko
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)

    for setup in mol_setups:
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
        if is_ok:
            with open(output_pdbqt_path, "w") as f:
                f.write(pdbqt_string)
            print(f"[+] Processed ligand PDBQT saved to: {output_pdbqt_path}")
            return
        else:
            print(f"[-] Error parameterizing {sdf_path.name}: {error_msg}")


def main():
    # Load frozen settings
    with open("config/docking_config.yaml", "r") as f:
        config = yaml.safe_load(f)

    processed_dir = Path("data/processed")
    processed_dir.mkdir(parents=True, exist_ok=True)

    # 1. Preprocess Receptors (.cif parsing)
    for target_name, raw_path in config["receptors"].items():
        raw_cif = Path(raw_path)
        if raw_cif.exists():
            clean_cif_path = processed_dir / f"{target_name}_clean.cif"
            clean_cif(raw_cif, clean_cif_path)
        else:
            print(f"[!] Warning: Raw file not found at {raw_path}")

    # 2. Preprocess Ligands (3D SDF -> PDBQT)
    ligand_dir = Path("data/raw/ligands")
    for ligand_name in config["ligands"]:
        sdf_file = ligand_dir / f"{ligand_name}.sdf"
        if sdf_file.exists():
            out_pdbqt = processed_dir / f"{ligand_name}.pdbqt"
            prepare_ligand_pdbqt(sdf_file, out_pdbqt)
        else:
            print(f"[!] Warning: Raw SDF file not found at {sdf_file}")


if __name__ == "__main__":
    main()