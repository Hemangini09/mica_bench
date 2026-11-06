#!/usr/bin/env python3
"""Run ModelAngelo on each B2 target when its MRC and FASTA exist."""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

from config import DATA, LOGS, MODELS, PANEL, RUN_PDB_IDS, ROOT

CONDA = Path("/home/hemangini/miniconda3/etc/profile.d/conda.sh")
TORCH_HOME = ROOT / "weights" / "torch"


def find_built_model(folder: Path):
    if not folder.exists():
        return None
    files = [
        p
        for p in folder.rglob("*")
        if p.suffix.lower() in {".cif", ".pdb"} and p.stat().st_size > 1000
    ]
    return max(files, key=lambda p: p.stat().st_size) if files else None


def inputs_ready(pdb: str) -> bool:
    info = PANEL[pdb]
    mrc = DATA / pdb / f"emd_{info['emdb']}.mrc"
    fasta = DATA / pdb / f"{pdb}.fasta"
    cif = DATA / pdb / f"{pdb}.cif"
    return all(p.exists() and p.stat().st_size > 100 for p in (mrc, fasta, cif))


def run_one(pdb: str) -> int:
    info = PANEL[pdb]
    out = MODELS / "modelangelo" / pdb
    if find_built_model(out):
        print(f"already built {pdb}", flush=True)
        return 0
    out.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    mrc = DATA / pdb / f"emd_{info['emdb']}.mrc"
    fasta = DATA / pdb / f"{pdb}.fasta"
    log = out / "run.log"
    script = f"""
set -e
source {CONDA}
conda activate model_angelo
export TORCH_HOME={TORCH_HOME}
export MPLBACKEND=Agg
model_angelo build -v {mrc} -f {fasta} -o {out} --device 0
"""
    print(f"ModelAngelo {pdb}", flush=True)
    with log.open("w") as handle:
        proc = subprocess.run(
            ["bash", "-lc", script],
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if proc.returncode != 0:
        print(f"FAIL {pdb} see {log}", flush=True)
    return proc.returncode


def main() -> int:
    wait = "--wait" in sys.argv
    failed = []
    remaining = list(RUN_PDB_IDS)
    while remaining:
        progressed = False
        still = []
        for pdb in remaining:
            if not inputs_ready(pdb):
                still.append(pdb)
                continue
            progressed = True
            if run_one(pdb) != 0:
                failed.append(pdb)
        remaining = still
        if not remaining:
            break
        if not wait:
            print("waiting on maps:", " ".join(remaining))
            break
        if not progressed:
            print("waiting for downloads:", " ".join(remaining), flush=True)
            time.sleep(60)
    print("failed", failed)
    return 1 if failed else 0


if __name__ == "__main__":
    os.environ.setdefault("TORCH_HOME", str(TORCH_HOME))
    sys.exit(main())
