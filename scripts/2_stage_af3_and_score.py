#!/usr/bin/env python3
"""Combine AF3 model_0 CIFs per PDB and score with 8_evaluate_metrics."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from config import DATA, MODELS, RUN_PDB_IDS

AA3 = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}


def split_cif(line: str) -> list[str]:
    parts, token, quote = [], [], None
    for char in line.strip():
        if quote:
            if char == quote:
                quote = None
            else:
                token.append(char)
        elif char in ("'", '"'):
            quote = char
        elif char.isspace():
            if token:
                parts.append("".join(token))
                token = []
        else:
            token.append(char)
    if token:
        parts.append("".join(token))
    return parts


def cif_ca_atoms(path: Path):
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    header, header_end = None, 0
    for i, line in enumerate(lines):
        if line.startswith("_atom_site."):
            header = [] if header is None else header
            header.append(line[len("_atom_site.") :].strip())
            header_end = i + 1
        elif header is not None:
            break
    if not header:
        return []
    index = {name: n for n, name in enumerate(header)}
    required = ("label_atom_id", "label_comp_id", "Cartn_x", "Cartn_y", "Cartn_z")
    if any(name not in index for name in required):
        return []
    atoms = []
    for line in lines[header_end:]:
        if not line.startswith("ATOM"):
            if atoms:
                break
            continue
        parts = split_cif(line)
        if len(parts) <= max(index[name] for name in required):
            continue
        if parts[index["label_atom_id"]] != "CA":
            continue
        atoms.append(
            (
                parts[index["label_comp_id"]],
                float(parts[index["Cartn_x"]]),
                float(parts[index["Cartn_y"]]),
                float(parts[index["Cartn_z"]]),
            )
        )
    return atoms


def write_pdb(atoms_by_chain, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    serial = 1
    lines = []
    for chain, atoms in atoms_by_chain:
        for i, (res, x, y, z) in enumerate(atoms, 1):
            res3 = res if len(res) == 3 else AA3.get(res, "UNK")
            lines.append(
                f"ATOM  {serial:5d}  CA  {res3:3s} {chain}{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n"
            )
            serial += 1
        lines.append("TER\n")
    lines.append("END\n")
    dest.write_text("".join(lines))


def main():
    chains_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for pdb in RUN_PDB_IDS:
        cifs = sorted((DATA / pdb / "AF3_chains").glob("*model_0.cif"))
        if not cifs:
            print("no AF3", pdb)
            continue
        packed = []
        for i, cif in enumerate(cifs):
            atoms = cif_ca_atoms(cif)
            packed.append((chains_ids[i], atoms))
            print(pdb, cif.name, "CA", len(atoms))
        write_pdb(packed, MODELS / "af3_alone" / pdb / "af3_combined.pdb")
    subprocess.check_call([sys.executable, str(Path(__file__).parent / "8_evaluate_metrics.py")])


if __name__ == "__main__":
    main()
