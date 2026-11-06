#!/usr/bin/env python3
"""Stage MICA inputs, process AF3 (CPU/Merizo), Phenix-dock into maps.

Runs entirely on CPU so it can overlap the Cryo2Struct GPU job.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

from config import DATA, FAST_FIRST, LOGS, MICA_SRC, PANEL, PHENIX_ACT, WORK

UA = {"User-Agent": "mica-bench/1.0"}
CONDA = Path("/home/hemangini/miniconda3/etc/profile.d/conda.sh")


def fetch_contour(emdb: str) -> float | None:
    url = f"https://www.ebi.ac.uk/emdb/api/entry/{emdb}"
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read().decode())
        return float(data["map"]["contour_list"]["contour"][0]["level"])
    except Exception as exc:
        print("contour fail", emdb, exc, flush=True)
        return None


def percentile_contour(mrc_path: Path) -> float:
    import mrcfile
    import numpy as np

    with mrcfile.open(mrc_path, permissive=True) as mrc:
        data = np.array(mrc.data, dtype=np.float32)
    nz = data[data > 0]
    if nz.size == 0:
        return 0.0
    return float(np.percentile(nz, 90))


def wait_for(path: Path, label: str, timeout_s: int = 7200) -> None:
    t0 = time.time()
    while not path.exists():
        if time.time() - t0 > timeout_s:
            raise SystemExit(f"timeout waiting for {label}: {path}")
        print("waiting", label, flush=True)
        time.sleep(20)


def mica_python() -> list[str]:
    env_py = Path("/home/hemangini/miniconda3/envs/MICA/bin/python")
    if env_py.exists():
        return [str(env_py)]
    return ["bash", "-lc", f"source {CONDA} && conda activate MICA && python"]


def stage_target(pdb: str) -> Path:
    info = PANEL[pdb]
    dest = WORK / pdb
    af3 = dest / "AF3_results"
    af3.mkdir(parents=True, exist_ok=True)
    fasta_src = DATA / pdb / f"{pdb}.fasta"
    shutil.copy2(fasta_src, dest / f"{pdb}.fasta")
    src_map = DATA / pdb / f"emd_{info['emdb']}.mrc"
    map_dest = dest / f"emd_{info['emdb']}.mrc"
    if not map_dest.exists():
        os.symlink(src_map, map_dest)
    jobs = sorted((DATA / pdb / "AF3").glob("*_u*"))
    fasta_ids = []
    header = None
    for line in fasta_src.read_text().splitlines():
        if line.startswith(">"):
            fasta_ids.append(line[1:].split("|")[0].strip())
    if not jobs:
        raise FileNotFoundError(f"no AF3 jobs for {pdb}")
    for i, job in enumerate(jobs):
        rec = fasta_ids[i] if i < len(fasta_ids) else f"{pdb}_{i+1}"
        cif = next(job.glob("*_model_0.cif"), None)
        if cif is None:
            print("no model_0", job, flush=True)
            continue
        out_dir = af3 / rec
        out_dir.mkdir(exist_ok=True)
        named = out_dir / f"fold_{rec.lower()}_model_0.cif"
        if named.exists() or named.is_symlink():
            continue
        os.symlink(cif.resolve(), named)
    return dest


def already_docked(dest: Path) -> bool:
    hits = list((dest / "AF3_docked_models").glob("*.pdb"))
    return any(p.stat().st_size > 500 for p in hits)


def run_process_af3(dest: Path, pdb: str) -> None:
    af3 = dest / "AF3_results"
    marker = dest / "AF3_domains"
    if marker.exists() and any(marker.glob("*.pdb")):
        print("af3 processed", pdb, flush=True)
        return
    fasta = dest / f"{pdb}.fasta"
    log = LOGS / f"mica_process_af3_{pdb}.log"
    LOGS.mkdir(parents=True, exist_ok=True)
    cmd = (
        f"source {CONDA} && conda activate MICA && "
        f"cd {MICA_SRC} && python utils/process_AF3_results.py "
        f"-f {fasta} -a {af3} -d cpu"
    )
    print("process_AF3", pdb, flush=True)
    with log.open("w") as handle:
        proc = subprocess.run(["bash", "-c", cmd], stdout=handle, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"process_AF3 failed {pdb} see {log}")


def run_dock(dest: Path, pdb: str, contour: float, resolution: float) -> None:
    if already_docked(dest):
        print("already docked", pdb, flush=True)
        return
    info = PANEL[pdb]
    log = LOGS / f"mica_dock_{pdb}.log"
    map_path = dest / f"emd_{info['emdb']}.mrc"
    fasta = dest / f"{pdb}.fasta"
    af3 = dest / "AF3_results"
    domains = dest / "AF3_domains"
    dock_type = "domainwise" if any(domains.glob("*.pdb")) else "chainwise"
    cmd = (
        f"source {CONDA} && conda activate MICA && "
        f"cd {MICA_SRC} && python utils/dock_in_map.py "
        f"-m {map_path} -c {contour} -r {resolution} "
        f"-f {fasta} -a {af3} --dock_type {dock_type} -x {PHENIX_ACT}"
    )
    print("dock", pdb, "contour", contour, "res", resolution, flush=True)
    with log.open("w") as handle:
        proc = subprocess.run(["bash", "-c", cmd], stdout=handle, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        print("FAIL dock", pdb, "see", log, flush=True)
    else:
        print("ok dock", pdb, flush=True)


def main() -> int:
    wait_for(Path("/home/hemangini/miniconda3/envs/MICA/bin/python"), "MICA env python", timeout_s=3600)
    wait_for(PHENIX_ACT, "phenix_env.sh")
    contours_path = WORK / "contours.json"
    WORK.mkdir(parents=True, exist_ok=True)
    contours: dict[str, float] = {}
    if contours_path.exists():
        contours = json.loads(contours_path.read_text())
    failed = []
    for pdb in FAST_FIRST:
        try:
            dest = stage_target(pdb)
            emdb = PANEL[pdb]["emdb"]
            if pdb not in contours:
                c = fetch_contour(emdb)
                if c is None:
                    c = percentile_contour(DATA / pdb / f"emd_{emdb}.mrc")
                    print("percentile contour", pdb, c, flush=True)
                contours[pdb] = c
                contours_path.write_text(json.dumps(contours, indent=2))
            run_process_af3(dest, pdb)
            run_dock(dest, pdb, contours[pdb], float(PANEL[pdb]["resolution"]))
        except Exception as exc:
            print("FAIL", pdb, exc, flush=True)
            failed.append(pdb)
    print("dock failed", failed, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
