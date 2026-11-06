#!/usr/bin/env python3
"""Fast-wave MICA: AF3_struct_free only (no Phenix dock / refine)."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from config import DATA, FAST_FIRST, LOGS, MICA_SRC, MODELS, PANEL, PULCHRA, WORK

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


def already(pdb: str) -> bool:
    for folder in (MODELS / "mica" / pdb, MODELS / "mica_struct_free" / pdb):
        if not folder.exists():
            continue
        # 0-byte all-atom files mean MICA finished with no aligned fragments.
        if any(
            p.name.endswith("_all_atom_model.pdb") and "backup" not in p.parts
            for p in folder.rglob("*.pdb")
        ):
            return True
        hits = [
            p
            for p in folder.rglob("*.pdb")
            if "(init)" not in p.name
            and "backup" not in p.parts
            and p.stat().st_size > 1000
        ]
        if hits:
            return True
    return False


def stage(pdb: str) -> Path:
    dest = WORK / pdb
    af3 = dest / "AF3_results"
    af3.mkdir(parents=True, exist_ok=True)
    fasta_src = DATA / pdb / f"{pdb}.fasta"
    shutil.copy2(fasta_src, dest / f"{pdb}.fasta")
    emdb = PANEL[pdb]["emdb"]
    src_map = DATA / pdb / f"emd_{emdb}.mrc"
    map_dest = dest / f"emd_{emdb}.mrc"
    if not map_dest.exists():
        os.symlink(src_map, map_dest)
    jobs = sorted((DATA / pdb / "AF3").glob("*_u*"))
    fasta_ids = []
    for line in fasta_src.read_text().splitlines():
        if line.startswith(">"):
            fasta_ids.append(line[1:].split("|")[0].strip())
    for i, job in enumerate(jobs):
        rec = fasta_ids[i] if i < len(fasta_ids) else f"{pdb}_{i+1}"
        cif = next(job.glob("*_model_0.cif"), None)
        if cif is None:
            continue
        out_dir = af3 / rec
        out_dir.mkdir(exist_ok=True)
        named = out_dir / f"fold_{rec.lower()}_model_0.cif"
        if not named.exists():
            os.symlink(cif.resolve(), named)
    return dest


def run_one(pdb: str) -> int:
    dest = stage(pdb)
    info = PANEL[pdb]
    out_dir = MODELS / "mica_struct_free" / pdb
    out_dir.mkdir(parents=True, exist_ok=True)
    log = LOGS / f"mica_fast_{pdb}.log"
    map_path = dest / f"emd_{info['emdb']}.mrc"
    fasta = dest / f"{pdb}.fasta"
    af3 = dest / "AF3_results"
    cmd = (
        f"source {CONDA} && conda activate MICA && "
        f"cd {MICA_SRC} && python run.py --protocol AF3_struct_free "
        f"-m {map_path} -f {fasta} -a {af3} --device cuda -o {out_dir} "
        f"--run_pulchra --pulchra_path {PULCHRA} --quiet"
    )
    print("MICA struct_free", pdb, flush=True)
    with log.open("w") as handle:
        return subprocess.run(["bash", "-c", cmd], stdout=handle, stderr=subprocess.STDOUT, text=True).returncode


def main() -> int:
    LOGS.mkdir(parents=True, exist_ok=True)
    only = None
    if "--only" in sys.argv:
        idx = sys.argv.index("--only")
        if idx + 1 >= len(sys.argv):
            raise SystemExit("usage: 6_run_mica.py --only PDB")
        only = sys.argv[idx + 1].upper()
    targets = [only] if only else list(FAST_FIRST)
    failed = []
    for pdb in targets:
        if already(pdb):
            print("already", pdb, flush=True)
            continue
        while True:
            psout = subprocess.check_output(["ps", "-eo", "args="], text=True)
            other_gpu = False
            for line in psout.splitlines():
                if "6_run_mica" in line or "pgrep" in line:
                    continue
                if "atom_inference.py" in line or "amino_inference.py" in line:
                    other_gpu = True
                if "cryoatom build" in line or "cryo2struct.py" in line:
                    other_gpu = True
                if "protocol=temp_flex" in line or "protocol AF3_struct" in line:
                    other_gpu = True
            if gpu_mem_mib() <= 28000 and not other_gpu:
                break
            print("waiting GPU", gpu_mem_mib(), "busy" if other_gpu else "", flush=True)
            time.sleep(45)
        rc = run_one(pdb)
        if already(pdb):
            print("ok (model on disk)", pdb, flush=True)
        elif rc != 0:
            print("FAIL", pdb, "see", LOGS / f"mica_fast_{pdb}.log", flush=True)
            failed.append(pdb)
        else:
            print("ok", pdb, flush=True)
    print("mica fast failed", failed, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
