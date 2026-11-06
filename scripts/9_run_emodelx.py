#!/usr/bin/env python3
"""EModelX(+AF) on FAST_FIRST maps. Uses AF3 models as ranked_0.pdb templates."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from config import DATA, FAST_FIRST, LOGS, MODELS, PANEL, ROOT

CONDA = Path("/home/hemangini/miniconda3/etc/profile.d/conda.sh")
SRC = ROOT / "src" / "EModelX"
PULCHRA = SRC / "modules" / "pulchra304" / "src" / "pulchra"


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
            "cryoatom build",
            "run.py --protocol=temp",
        )
    )


def wait_gpu() -> None:
    while gpu_busy() or gpu_mem_mib() > 22000:
        print("waiting GPU", gpu_mem_mib(), flush=True)
        time.sleep(45)


def already(pdb: str) -> bool:
    folder = MODELS / "emodelx" / pdb
    if not folder.exists():
        return False
    return any(p.stat().st_size > 1000 for p in folder.rglob("*.pdb") if "(init)" not in p.name)


def cif_to_pdb(cif: Path, pdb_out: Path) -> None:
    text = cif.read_text(errors="ignore")
    atoms = []
    header, in_loop = [], False
    cols: dict[str, int] = {}
    for line in text.splitlines():
        if line.startswith("_atom_site."):
            header.append(line[len("_atom_site.") :].strip())
            in_loop = True
            continue
        if in_loop and (line.startswith("_") or line.startswith("loop_") or line.startswith("#")):
            if cols:
                break
            continue
        if in_loop and header and not line.startswith("_") and line.strip() and not line.startswith("loop_"):
            if not cols:
                cols = {n: i for i, n in enumerate(header)}
            parts = line.split()
            if len(parts) < len(header):
                continue
            try:
                group = parts[cols.get("group_PDB", 0)]
                if group not in ("ATOM", "HETATM"):
                    continue
                name = parts[cols["label_atom_id"]].strip('"')
                resn = parts[cols["label_comp_id"]]
                chain = parts[cols.get("auth_asym_id", cols.get("label_asym_id", 0))][0]
                seq = int(float(parts[cols.get("auth_seq_id", cols["label_seq_id"])]))
                x, y, z = (float(parts[cols[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z"))
                occ = float(parts[cols["occupancy"]]) if "occupancy" in cols else 1.0
                b = float(parts[cols["B_iso_or_equiv"]]) if "B_iso_or_equiv" in cols else 0.0
                elem = parts[cols["type_symbol"]] if "type_symbol" in cols else name[0]
            except (KeyError, ValueError, IndexError):
                continue
            atoms.append((name, resn, chain, seq, x, y, z, occ, b, elem))
    pdb_out.parent.mkdir(parents=True, exist_ok=True)
    with pdb_out.open("w") as handle:
        for i, (name, resn, chain, seq, x, y, z, occ, b, elem) in enumerate(atoms, 1):
            handle.write(
                f"ATOM  {i:5d} {name:>4s} {resn:3s} {chain}{seq:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}          {elem:>2s}\n"
            )
        handle.write("END\n")
    if pdb_out.stat().st_size < 200:
        raise RuntimeError(f"empty pdb from {cif}")


def stage_templates(pdb: str) -> Path:
    dest = ROOT / "work" / pdb / "emodelx_templates"
    dest.mkdir(parents=True, exist_ok=True)
    fasta = DATA / pdb / f"{pdb}.fasta"
    recs = []
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            recs.append(line[1:].split("|")[0].strip())
    jobs = sorted((DATA / pdb / "AF3").glob("*_u*"))
    for i, job in enumerate(jobs):
        cif = next(job.glob("*_model_0.cif"), None)
        if cif is None:
            continue
        rec = recs[i] if i < len(recs) else f"{pdb}_{i+1}"
        out = dest / rec / "ranked_0.pdb"
        if not out.exists() or out.stat().st_size < 200:
            cif_to_pdb(cif, out)
    return dest


def run_one(pdb: str) -> int:
    wait_gpu()
    info = PANEL[pdb]
    out = MODELS / "emodelx" / pdb
    out.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    templates = stage_templates(pdb)
    log = LOGS / f"emodelx_{pdb}.log"
    mrc = DATA / pdb / f"emd_{info['emdb']}.mrc"
    fasta = DATA / pdb / f"{pdb}.fasta"
    cmd = (
        f"source {CONDA} && conda activate EModelX && cd {SRC} && "
        f"python run.py --protocol=temp_flex --EM_map={mrc} --fasta={fasta} "
        f"--template_dir={templates} --output_dir={out} --run_pulchra "
        f"--pulchra_path={PULCHRA} --resolution={info['resolution']}"
    )
    print("EModelX", pdb, flush=True)
    with log.open("w") as handle:
        return subprocess.run(["bash", "-c", cmd], stdout=handle, stderr=subprocess.STDOUT, text=True).returncode


def main() -> int:
    failed = []
    for pdb in FAST_FIRST:
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
    print("emodelx failed", failed, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
