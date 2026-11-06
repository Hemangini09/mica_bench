#!/usr/bin/env python3
"""CryoAtom2 on FAST_FIRST maps. Waits until the GPU is free of other DL jobs."""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

from config import DATA, FAST_FIRST, LOGS, MODELS, PANEL, ROOT

CONDA = Path("/home/hemangini/miniconda3/etc/profile.d/conda.sh")


def gpu_mem_mib() -> int:
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
            text=True,
        )
        return int(out.strip().split()[0])
    except Exception:
        return 0


def gpu_busy() -> bool:
    psout = subprocess.check_output(["ps", "-eo", "args="], text=True)
    return any(
        tok in psout
        for tok in (
            "atom_inference.py",
            "amino_inference.py",
            "AF3_struct_free",
            "run.py --protocol=temp",
        )
    )


def wait_gpu() -> None:
    while gpu_busy() or gpu_mem_mib() > 22000:
        print("waiting GPU", gpu_mem_mib(), flush=True)
        time.sleep(45)


def already(pdb: str) -> bool:
    folder = MODELS / "cryoatom" / pdb
    if not folder.exists():
        return False
    return any(p.stat().st_size > 1000 for p in folder.rglob("*.pdb") if "dots" not in p.name)


def run_one(pdb: str) -> int:
    wait_gpu()
    info = PANEL[pdb]
    out = MODELS / "cryoatom" / pdb
    out.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    log = LOGS / f"cryoatom_{pdb}.log"
    mrc = DATA / pdb / f"emd_{info['emdb']}.mrc"
    fasta = DATA / pdb / f"{pdb}.fasta"
    cmd = (
        f"source {CONDA} && conda activate CryoAtom2 && export MPLBACKEND=Agg && "
        f"cryoatom build -v {mrc} -ps {fasta} -o {out} -d 0"
    )
    print("CryoAtom", pdb, flush=True)
    with log.open("w") as handle:
        return subprocess.run(["bash", "-c", cmd], stdout=handle, stderr=subprocess.STDOUT, text=True).returncode


def main() -> int:
    failed = []
    for pdb in FAST_FIRST:
        if pdb == "9TBI":
            print("skip 9TBI too large", flush=True)
            continue
        if already(pdb):
            print("already", pdb, flush=True)
            continue
        rc = run_one(pdb)
        if already(pdb):
            print("ok (model on disk)", pdb, flush=True)
        elif rc != 0:
            print("FAIL", pdb, flush=True)
            failed.append(pdb)
        else:
            print("ok", pdb, flush=True)
    print("cryoatom failed", failed, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
