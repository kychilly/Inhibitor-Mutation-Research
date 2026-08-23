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


def prepare_ligand_pdbqt(sdf_path, output_pdbqt_path, num_conformers=15):
    """
    Generates multiple 3D conformers, minimizes them via MMFF94 forcefield,
    extracts the lowest-energy conformer, and exports it to AutoDock PDBQT format.
    """
    supplier = Chem.SDMolSupplier(str(sdf_path))
    mol = next(supplier)
    if mol is None:
        raise ValueError(f"Failed to load ligand from {sdf_path}")

    # Add polar hydrogens required for forcefield calculations and docking
    mol = Chem.AddHs(mol, addCoords=True)

    # 1. Generate 15 distinct 3D rotamer conformers
    conf_ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=num_conformers,
        randomSeed=42,
        pruneRmsThresh=0.5
    )

    # 2. Minimize each conformer using MMFF94 forcefield
    minimized_energies = []
    for cid in conf_ids:
        AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=500)
        ff = AllChem.MMFFGetMoleculeForceField(
            mol, AllChem.MMFFGetMoleculeProperties(mol), confId=cid
        )
        if ff:
            minimized_energies.append((cid, ff.CalcEnergy()))

    if not minimized_energies:
        raise RuntimeError(f"MMFF94 energy minimization failed for {sdf_path.name}")

    # 3. Identify lowest-energy conformer ID
    best_conf_id = min(minimized_energies, key=lambda x: x[1])[0]
    best_energy = min(minimized_energies, key=lambda x: x[1])[1] # PE can be negative, Etotal = Ebonds + Eangles + Edihedrals + Evanderwalls + Eelectrostatic
    print(f"[+] {sdf_path.stem} | Best Conformer ID: {best_conf_id} | Potential Energy(yes this can be negative): {best_energy:.2f} kcal/mol")

    # 4. Isolate ONLY the lowest-energy conformer into a single-conformer Mol
    best_conf = mol.GetConformer(best_conf_id)
    single_conf_mol = Chem.Mol(mol)
    single_conf_mol.RemoveAllConformers()
    single_conf_mol.AddConformer(best_conf, assignId=True)

    # 5. Convert to PDBQT via Meeko
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(single_conf_mol)

    for setup in mol_setups:
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
        if is_ok:
            with open(output_pdbqt_path, "w") as f:
                f.write(pdbqt_string)
            print(f"[+] Processed & saved lowest-energy ligand PDBQT: {output_pdbqt_path}")
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