#!/usr/bin/env python3
"""Run phenix.map_to_model on each B2 target (CPU; EMProt replacement)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from config import DATA, FAST_FIRST, LOGS, MODELS, PANEL, PHENIX_ACT

CONDA = Path("/home/hemangini/miniconda3/etc/profile.d/conda.sh")


def already(pdb: str) -> bool:
    folder = MODELS / "phenix_map_to_model" / pdb
    if not folder.exists():
        return False
    hit = folder / "map_to_model.pdb"
    return hit.exists() and hit.stat().st_size > 1000


def run_one(pdb: str) -> int:
    if already(pdb):
        print("already", pdb, flush=True)
        return 0
    info = PANEL[pdb]
    out = MODELS / "phenix_map_to_model" / pdb
    out.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    mrc = DATA / pdb / f"emd_{info['emdb']}.mrc"
    fasta = DATA / pdb / f"{pdb}.fasta"
    log = LOGS / f"map_to_model_{pdb}.log"
    pdb_out = out / "map_to_model.pdb"
    cmd = (
        f"source {PHENIX_ACT} && cd {out} && "
        f"phenix.map_to_model {mrc} {fasta} "
        f"resolution={info['resolution']} nproc=16 "
        f"pdb_out={pdb_out.name} scattering_table=electron"
    )
    print("map_to_model", pdb, flush=True)
    with log.open("w") as handle:
        proc = subprocess.run(["bash", "-c", cmd], stdout=handle, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0 or not pdb_out.exists() or pdb_out.stat().st_size < 1000:
        hits = [p for p in out.rglob("*.pdb") if p.stat().st_size > 1000]
        if hits:
            print("ok map_to_model (alt pdb)", pdb, hits[0].name, flush=True)
            return 0
        print("FAIL map_to_model", pdb, "see", log, flush=True)
        return proc.returncode or 1
    print("ok map_to_model", pdb, flush=True)
    return 0


def main() -> int:
    failed = []
    for pdb in FAST_FIRST:
        try:
            if run_one(pdb) != 0:
                failed.append(pdb)
        except Exception as exc:
            print("FAIL", pdb, exc, flush=True)
            failed.append(pdb)
    print("map_to_model failed", failed, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
