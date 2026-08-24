#!/usr/bin/env python3
"""
organize_plip_files.py

Sorts intermediate PLIP output files (e.g., plipfixed.<receptor>_<drug>_complex_<hash>.pdb)
located directly under results/plip_results/ into their corresponding target subdirectories:
results/plip_results/{receptor}/{drug}/
"""

import re
import shutil
from pathlib import Path

PLIP_RESULTS_DIR = Path("results/plip_results")
RECEPTORS = ["wildtype", "t790m"]
LIGANDS = ["afatinib", "dacomitinib", "erlotinib", "gefitinib", "osimertinib"]


def organize_files():
    if not PLIP_RESULTS_DIR.exists():
        print(f"[-] Directory not found: {PLIP_RESULTS_DIR}")
        return

    # Find all files directly under results/plip_results/
    loose_files = [f for f in PLIP_RESULTS_DIR.iterdir() if f.is_file()]

    if not loose_files:
        print("[+] No loose files found in results/plip_results to organize.")
        return

    moved_count = 0

    for file_path in loose_files:
        filename = file_path.name

        # Skip summary text files or other non-pdb artifacts
        if filename == "hbond_loss_summary.txt":
            continue

        # Match receptor and ligand from pattern: plipfixed.<receptor>_<drug>_complex...
        match = re.search(r"plipfixed\.([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_complex", filename)
        if match:
            receptor, ligand = match.group(1).lower(), match.group(2).lower()

            if receptor in RECEPTORS and ligand in LIGANDS:
                target_dir = PLIP_RESULTS_DIR / receptor / ligand
                target_dir.mkdir(parents=True, exist_ok=True)

                dest_path = target_dir / filename
                shutil.move(str(file_path), str(dest_path))
                print(f"  [Moved] {filename} -> {receptor}/{ligand}/")
                moved_count += 1
            else:
                print(f"  [Skipped] Unknown receptor/ligand in filename: {filename}")
        else:
            print(f"  [Skipped] Pattern mismatch: {filename}")

    print(f"\n[+] File organization complete. Moved {moved_count} file(s).")


if __name__ == "__main__":
    organize_files()