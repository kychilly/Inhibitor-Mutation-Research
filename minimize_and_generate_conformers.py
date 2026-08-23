# might be useless file, remove later if necessary

from rdkit import Chem
from rdkit.Chem import AllChem


def minimize_and_generate_conformers(sdf_path, num_conformers=15):
   # Generates 15 distinct 3D conformers for a ligand and applies
   # MMFF94 forcefield energy minimization to find the most stable pose.
    suppl = Chem.SDMolSupplier(str(sdf_path))
    mol = next(suppl)
    mol = Chem.AddHs(mol)  # Add explicit hydrogens for forcefield calculation

    # 1. Embed 10 distinct 3D rotamer conformers
    conformer_ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=num_conformers,
        randomSeed=42,
        pruneRmsThresh=0.5
    )

    # 2. Minimize energy for each conformer using MMFF94
    minimized_energies = []
    for conf_id in conformer_ids:
        # Returns (check_code, energy_value)
        results = AllChem.MMFFOptimizeMoleculeConformer(mol, confId=conf_id, maxIters=500)

        # Calculate forcefield potential energy
        ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id)
        if ff:
            minimized_energies.append((conf_id, ff.CalcEnergy()))

    # 3. Select the conformer with the lowest potential energy (global minimum pose)
    best_conf_id = min(minimized_energies, key=lambda x: x[1])[0]

    return mol, best_conf_id