import os
import glob
from Bio.PDB import MMCIFParser, PDBIO, Select
from openbabel import openbabel as ob

# This method essentially only serves to put the receptors into the processed data folder
# while also turning them into pdbqt files instead of cif for the docking_runner method to read/run through

RAW_RECEPTORS_DIR = "data/raw/receptors"
PROCESSED_RECEPTORS_DIR = "data/processed/receptors"


class NonWaterSelect(Select):
    """Filter out water molecules and non-protein heteroatoms from the CIF file."""

    def accept_residue(self, residue):
        # Reject water molecules (HOH / WAT)
        return residue.get_resname() not in ["HOH", "WAT"]


def cif_to_pdb(cif_path, pdb_path):
    """Converts a raw CIF file into a cleaned PDB structure."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("receptor", cif_path)

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path, select=NonWaterSelect())
    print(f"  [Converted] {cif_path} -> {pdb_path} (Stripped waters)")


def pdb_to_pdbqt(pdb_path, pdbqt_path):
    """
    Adds polar hydrogens, assigns Gasteiger charges,
    and converts PDB to PDBQT using OpenBabel.
    """
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdbqt")

    mol = ob.OBMol()
    if not obConversion.ReadFile(mol, pdb_path):
        raise RuntimeError(f"OpenBabel failed to read {pdb_path}")

    # Add hydrogens suited for pH 7.4 (polar hydrogens essential for Vina)
    mol.AddHydrogens(False, True, 7.4)

    if not obConversion.WriteFile(mol, pdbqt_path):
        raise RuntimeError(f"OpenBabel failed to write {pdbqt_path}")

    print(f"  [Prepared] {pdb_path} -> {pdbqt_path} (Added polar H + Gasteiger charges)")


def process_all_receptors():
    os.makedirs(PROCESSED_RECEPTORS_DIR, exist_ok=True)

    cif_files = glob.glob(os.path.join(RAW_RECEPTORS_DIR, "*.cif"))
    if not cif_files:
        print(f"No CIF files found in {RAW_RECEPTORS_DIR}")
        return

    print("=== Processing Raw Receptors (CIF -> PDBQT) ===")
    for cif_path in cif_files:
        base_name = os.path.splitext(os.path.basename(cif_path))[0]

        # Map IDs to standardized processed filenames expected by docking_runner
        if base_name == "1M17":
            out_filename = "1M17_prepared.pdbqt"
        elif base_name == "2JIT":
            out_filename = "2JIT_prepared.pdbqt"
        else:
            out_filename = f"{base_name}_prepared.pdbqt"

        temp_pdb = os.path.join(PROCESSED_RECEPTORS_DIR, f"{base_name}_temp.pdb")
        final_pdbqt = os.path.join(PROCESSED_RECEPTORS_DIR, out_filename)

        # 1. Strip waters & convert CIF to standard PDB
        cif_to_pdb(cif_path, temp_pdb)

        # 2. Add polar H, assign charges, and output PDBQT
        pdb_to_pdbqt(temp_pdb, final_pdbqt)

        # Clean up temporary PDB file
        if os.path.exists(temp_pdb):
            os.remove(temp_pdb)

    print("\nReceptor preprocessing complete! Files ready in 'data/processed/receptors/'.")


if __name__ == "__main__":
    process_all_receptors()