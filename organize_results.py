# Feel free to incorporate this into docking_runner as it might(it is) more efficient
# However, you can also just run py organize_results.py after docking_runner has finished running to organize the files

import os
import re
import shutil
from pathlib import Path


def organize_docking_results(base_dir="data/processed/docking_results"):
    base_path = Path(base_dir)

    if not base_path.exists():
        print(f"[-] Directory not found: {base_path}")
        return

    # Regular expression to parse filenames: {receptor}_{ligand}_seed{num}.{ext}
    file_pattern = re.compile(r"^(wildtype|t790m)_([a-zA-Z0-0_-]+)_seed\d+\.(log|pdbqt)$")

    moved_count = 0

    for file_path in base_path.iterdir():
        if file_path.is_file():
            match = file_pattern.match(file_path.name)
            if match:
                receptor, ligand, extension = match.groups()

                # Determine target subfolder ('logs' or 'pdbqt')
                file_type_folder = "logs" if extension == "log" else "pdbqt"

                # Construct nested target directory path
                target_dir = base_path / receptor / ligand / file_type_folder
                target_dir.mkdir(parents=True, exist_ok=True)

                # Move file to target directory
                shutil.move(str(file_path), str(target_dir / file_path.name))
                moved_count += 1

    print(f"[+] Successfully organized {moved_count} docking result files.")


if __name__ == "__main__":
    organize_docking_results()