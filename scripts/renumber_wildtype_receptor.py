#!/usr/bin/env python3
"""
renumber_wildtype_receptor.py

1M17 (wildtype) uses author-numbering offset by -24 relative to the
canonical/UniProt EGFR numbering used by the T790M mutant structure --
confirmed by comparing PLIP contact residues across conditions (e.g.
wildtype MET769 and t790m's MET793 are the same physical hinge residue,
wildtype CYS773 / t790m CYS797 are the same covalent-warhead cysteine,
wildtype THR766 / t790m THR790 are the same gatekeeper position -- all a
clean +24 shift).

This script renumbers every ATOM/HETATM residue sequence number (PDB
columns 23-26) in the wildtype receptor PDB by +24, so both receptors
report contacts in the same numbering frame before going through PLIP.

The original file is backed up alongside the corrected one (suffix
".preoffset.pdb") rather than silently discarded, in case anything needs
to be re-derived from the raw author numbering later.

Usage:
    python scripts/renumber_wildtype_receptor.py
"""

from pathlib import Path

RECEPTOR_DIR = Path("data/processed/receptors")
WILDTYPE_PDB = RECEPTOR_DIR / "wildtype_clean.pdb"
OFFSET = 24


def renumber_residues(in_path: Path, out_path: Path, offset: int):
    """
    Adds `offset` to the residue sequence number (PDB columns 23-26,
    0-indexed [22:26]) on every ATOM/HETATM line. All other fields --
    including atom serials, names, coordinates, occupancy, and element --
    are left untouched.
    """
    lines_out = []
    n_atoms_shifted = 0
    seen_resnums = set()

    with open(in_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 26:
                try:
                    old_resnum = int(line[22:26])
                except ValueError:
                    # Non-numeric or blank resSeq field -- leave line untouched
                    lines_out.append(line)
                    continue

                new_resnum = old_resnum + offset
                seen_resnums.add(old_resnum)

                if new_resnum > 9999:
                    raise ValueError(
                        f"Renumbered residue {new_resnum} exceeds the 4-character "
                        f"PDB resSeq field width -- offset of {offset} is too large "
                        f"for this structure's numbering range."
                    )

                new_line = line[:22] + f"{new_resnum:>4}" + line[26:]
                lines_out.append(new_line)
                n_atoms_shifted += 1
            else:
                lines_out.append(line)

    with open(out_path, "w", encoding="utf-8") as f:
        f.writelines(lines_out)

    return n_atoms_shifted, (min(seen_resnums), max(seen_resnums)) if seen_resnums else (None, None)


def main():
    if not WILDTYPE_PDB.exists():
        print(f"[!] FATAL: wildtype receptor not found at {WILDTYPE_PDB.resolve()}")
        print("    Nothing to renumber -- check the path or rerun the .cif -> .pdb conversion first.")
        return

    backup_path = WILDTYPE_PDB.with_suffix(".preoffset.pdb")

    if backup_path.exists():
        print(f"[i] Backup already exists at {backup_path.resolve()} -- not overwriting it.")
        print("    Delete it manually first if you want to re-run this from the original numbering.")
        return

    # Preserve the pre-offset original before we overwrite wildtype_clean.pdb
    WILDTYPE_PDB.rename(backup_path)
    print(f"[i] Backed up original (author-numbered) file to {backup_path.resolve()}")

    n_shifted, (old_min, old_max) = renumber_residues(backup_path, WILDTYPE_PDB, OFFSET)

    print(f"[+] Renumbered {n_shifted} atom lines in {WILDTYPE_PDB.resolve()}")
    if old_min is not None:
        print(f"    Residue range shifted: {old_min}-{old_max}  ->  {old_min + OFFSET}-{old_max + OFFSET}")
    print(f"[+] Done. {WILDTYPE_PDB.name} now uses canonical EGFR numbering "
          f"(matches t790m_clean.pdb's frame).")


if __name__ == "__main__":
    main()