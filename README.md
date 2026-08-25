Markdown
## License

This project is open-source and available under the [MIT License](LICENSE).

Computational modeling and structural profiling pipeline evaluating 1st, 2nd, and 3rd generation EGFR tyrosine kinase inhibitors (TKIs) against Wildtype (WT) and T790M gatekeeper mutant receptors.

## Data Registration Log

**Access Date:** 2026-08-21

### Target Structures & Ligand Library Overview

| Structural Target / Dataset        | Accession / Source    | Resolution / Format | Local Directory                     |
|:-----------------------------------|:----------------------| :--- |:------------------------------------|
| **EGFR Kinase Domain (WT)**        | RCSB PDB: `1M17`      | 2.60 Å (PDB / PDBQT) | `data/raw/receptors/1M17.cif`       |
| **EGFR T790M Mutant**              | RCSB PDB: `224B`      | 2.80 Å (PDB / PDBQT) | `data/raw/receptors/2JIT.cif/`      |
| **Gefitinib: 1st Gen Inhibitor**   | PubChem (Gefitinib)   | SDF / PDBQT | `data/raw/ligands/gefitinib.sdf/`   |
| **Erlotinib: 1st Gen Inhibitor**   | PubChem (Erlotinib)   | SDF / PDBQT | `data/raw/ligands/erlotinib.sdf`    |
| **Afatinib: 2nd Gen Inhibitor**    | PubChem (Afatinib)    | SDF / PDBQT | `data/raw/ligands/afatinib.sdf/`    |
| **Dacomitinib: 2nd Gen Inhibitor** | PubChem (Dacomitinib) | SDF / PDBQT | `data/raw/ligands/acomitinib.sdf/`  |
| **Osimertinib: 3rd Gen Inhibitor** | PubChem (Osimertinib) | SDF / PDBQT | `data/raw/ligands/osimertinib.sdf/` |

## For replication, pull and follow:
Execution Order

1. python scripts/docking_runner.py          - Execute docking simulations

2. python scripts/plip_analyzer.py           - Analyze interaction fingerprints

3. python scripts/build_results_table.py     - Build summary data tables

4. python scripts/plot_results.py            - Generate binding affinity plots

5. python scripts/generate_pymol_figures.py  - Render 3D PyMOL structural models

6. python scripts/statistical_validation.py  - Perform 2D diagram generation & IC50 validation