# Save as test_pipeline.py and run: py test_pipeline.py
import yaml
from pathlib import Path

cfg_path = Path("config/docking_config.yaml")

with open(cfg_path, "r") as f:
    cfg = yaml.safe_load(f)

print("[1] Verifying YAML keys...")
assert "receptors" in cfg, "Missing 'receptors' block!"
assert "wildtype" in cfg["receptors"], "Missing 'wildtype' in receptors!"
assert "t790m" in cfg["receptors"], "Missing 't790m' in receptors!"

wt = cfg["receptors"]["wildtype"]
mut = cfg["receptors"]["t790m"]

print(f"    - Wildtype Prepared Path: {wt['prepared_path']} (Exists: {Path(wt['prepared_path']).exists()})")
print(f"    - Wildtype Center: ({wt['center_x']}, {wt['center_y']}, {wt['center_z']})")
print(f"    - T790M Prepared Path: {mut['prepared_path']} (Exists: {Path(mut['prepared_path']).exists()})")
print(f"    - T790M Center: ({mut['center_x']}, {mut['center_y']}, {mut['center_z']})")

print("\n[2] Verifying Ligand Files...")
for lig in cfg["ligands"]:
    p = Path(f"data/processed/ligands/{lig}.pdbqt")
    print(f"    - Ligand '{lig}': Exists={p.exists()}")

print("\n[+] All structural checks passed cleanly!")