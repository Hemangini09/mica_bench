#!/usr/bin/env python3
"""Download EMDB maps, RCSB FASTA, and deposited CIFs for RUN_PDB_IDS."""

from __future__ import annotations

import gzip
import shutil
import sys
import urllib.request
from pathlib import Path

from config import DATA, PANEL, RUN_PDB_IDS

UA = {"User-Agent": "mica-bench/1.0 (research; local benchmark)"}


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        print(f"skip {dest.name}", flush=True)
        return dest
    print(f"download {dest.name}", flush=True)
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=300) as src, dest.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    return dest


def unique_fasta(text: str) -> str:
    records = []
    header, seq = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(seq)))
            header, seq = line, []
        else:
            seq.append(line.strip())
    if header:
        records.append((header, "".join(seq)))
    seen, out = set(), []
    for header, seq in records:
        if seq in seen:
            continue
        seen.add(seq)
        out.append(header + "\n" + seq + "\n")
    return "".join(out)


def download_pdb(pdb: str) -> None:
    info = PANEL[pdb]
    folder = DATA / pdb
    folder.mkdir(parents=True, exist_ok=True)
    number = info["emdb"]
    gz = fetch(
        f"https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-{number}/map/emd_{number}.map.gz",
        folder / f"emd_{number}.map.gz",
    )
    mrc = folder / f"emd_{number}.mrc"
    if not mrc.exists() or mrc.stat().st_size < 1000:
        print(f"decompress {mrc.name}", flush=True)
        with gzip.open(gz, "rb") as src, mrc.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    fasta = folder / f"{pdb}.fasta"
    if not fasta.exists() or fasta.stat().st_size < 10:
        print(f"fasta {pdb}", flush=True)
        req = urllib.request.Request(
            f"https://www.rcsb.org/fasta/entry/{pdb}/download", headers=UA
        )
        raw = urllib.request.urlopen(req, timeout=120).read().decode()
        fasta.write_text(unique_fasta(raw))
    fetch(f"https://files.rcsb.org/download/{pdb}.cif", folder / f"{pdb}.cif")


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    failed = []
    for pdb in RUN_PDB_IDS:
        try:
            download_pdb(pdb)
            print(f"ready {pdb}", flush=True)
        except Exception as exc:
            failed.append((pdb, str(exc)))
            print(f"FAIL {pdb}: {exc}", flush=True)
    print(f"done {len(RUN_PDB_IDS) - len(failed)}/{len(RUN_PDB_IDS)}")
    for pdb, err in failed:
        print(f"  {pdb}: {err}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
