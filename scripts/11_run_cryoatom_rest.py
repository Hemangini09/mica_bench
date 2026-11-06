#!/usr/bin/env python3
"""CryoAtom on the small maps it has not built yet. Skips 9TBI."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from config import DATA, LOGS, MODELS, PANEL

# Smallest unfinished maps. 9XON/9W22/9WSD already have models. 9TBI is excluded.
ORDER = ("9N5N", "33AW", "9UIT", "9S3M")


def has_model(pdb: str) -> bool:
    folder = MODELS / "cryoatom" / pdb
    if not folder.exists():
        return False
    return any(
        p.stat().st_size > 1000 and "points" not in p.name
        for p in folder.rglob("*.cif")
    )


def run_one(pdb: str) -> int:
    info = PANEL[pdb]
    out = MODELS / "cryoatom" / pdb
    out.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    mrc = DATA / pdb / f"emd_{info['emdb']}.mrc"
    fasta = DATA / pdb / f"{pdb}.fasta"
    log = LOGS / f"cryoatom_{pdb}.log"
    cmd = f"cryoatom build -v {mrc} -ps {fasta} -o {out} -d 0"
    print("CryoAtom", pdb, flush=True)
    with log.open("w") as handle:
        return subprocess.run(
            ["bash", "-c", cmd], stdout=handle, stderr=subprocess.STDOUT
        ).returncode


def main() -> int:
    failed = []
    for pdb in ORDER:
        if has_model(pdb):
            print("already", pdb, flush=True)
            continue
        rc = run_one(pdb)
        ok = has_model(pdb)
        print("rc", pdb, rc, "model", ok, flush=True)
        if not ok:
            failed.append(pdb)
            break
    print("cryoatom failed", failed, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
