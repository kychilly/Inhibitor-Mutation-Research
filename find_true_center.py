# Save as find_true_center.py and run: py find_true_center.py
import numpy as np

def get_protein_geometric_center(pdbqt_path):
    coords = []
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except (ValueError, IndexError):
                    continue
    if not coords:
        print(f"[-] No ATOM records found in {pdbqt_path}")
        return None
    center = np.mean(coords, axis=0)
    return [round(float(c), 3) for c in center]

print("2JIT (T790M) Center: ", get_protein_geometric_center("data/processed/receptors/2JIT_prepared.pdbqt"))
print("1M17 (WT) Center:    ", get_protein_geometric_center("data/processed/receptors/1M17_prepared.pdbqt"))