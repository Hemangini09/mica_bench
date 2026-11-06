#!/usr/bin/env python3
"""Resample maps to 1 Å (no ChimeraX) and run Cryo2Struct on each B2 target."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import zoom

from config import DATA, FAST_FIRST, LOGS, MODELS, PANEL, ROOT

C2S = ROOT / "src" / "Cryo2Struct"
INPUT = C2S / "input"
CONDA = Path("/home/hemangini/miniconda3/etc/profile.d/conda.sh")


def resample_normalize(src_mrc: Path, dest_norm: Path) -> None:
    import mrcfile

    with mrcfile.open(src_mrc, permissive=True) as mrc:
        data = np.array(mrc.data, dtype=np.float32)
        voxel = np.array(mrc.voxel_size.tolist(), dtype=np.float64)
        origin = np.array(mrc.header.origin.tolist(), dtype=np.float64)
    # voxel_size may be a scalar-like
    if np.ndim(voxel) == 0:
        spacing = float(voxel)
        factors = (spacing, spacing, spacing)
    else:
        factors = tuple(float(v) / 1.0 for v in voxel[::-1])  # z,y,x vs x,y,z
        # mrcfile voxel_size is x,y,z; data is typically z,y,x
        vx, vy, vz = (float(voxel[0]), float(voxel[1]), float(voxel[2]))
        factors = (vz, vy, vx)
    if min(abs(f) for f in factors) < 1e-6:
        raise RuntimeError(f"bad voxel size {voxel} for {src_mrc}")
    resampled = zoom(data, factors, order=1).astype(np.float32)
    nz = resampled[np.nonzero(resampled)]
    if nz.size:
        resampled = resampled / np.percentile(nz, 95)
    resampled[resampled < 0] = 0
    resampled[resampled > 1] = 1
    dest_norm.parent.mkdir(parents=True, exist_ok=True)
    with mrcfile.new(dest_norm, overwrite=True) as out:
        out.set_data(resampled)
        out.voxel_size = 1
        try:
            ox, oy, oz = (float(origin[0]), float(origin[1]), float(origin[2]))
            out.header.origin.x = ox
            out.header.origin.y = oy
            out.header.origin.z = oz
        except Exception:
            pass


def stage_target(pdb: str) -> str:
    info = PANEL[pdb]
    emdb = info["emdb"]
    folder = INPUT / emdb
    folder.mkdir(parents=True, exist_ok=True)
    src_mrc = DATA / pdb / f"emd_{emdb}.mrc"
    fasta_src = DATA / pdb / f"{pdb}.fasta"
    shutil.copy2(fasta_src, folder / f"{pdb}.fasta")
    shutil.copy2(src_mrc, folder / f"emd_{emdb}.map")
    norm = folder / "emd_normalized_map.mrc"
    if not norm.exists() or norm.stat().st_size < 1000:
        print("resample", pdb, flush=True)
        resample_normalize(src_mrc, norm)
    return emdb


def already_built(emdb: str) -> bool:
    folder = INPUT / emdb
    hits = list(folder.glob("*cryo2struct*.pdb")) + list(folder.glob("*.pdb"))
    return any(p.stat().st_size > 1000 for p in hits)


def run_one(pdb: str, emdb: str) -> int:
    out_dir = MODELS / "cryo2struct" / pdb
    if already_built(emdb) or (out_dir.exists() and any(out_dir.rglob("*.pdb"))):
        print("already", pdb, flush=True)
        return 0
    LOGS.mkdir(parents=True, exist_ok=True)
    log = LOGS / f"cryo2struct_{pdb}.log"
    yml = C2S / "config" / "arguments.yml"
    yml.write_text(
        "\n".join(
            [
                f"density_map_name: '{emdb}'",
                f"input_data_dir: '{INPUT}'",
                f"amino_checkpoint: '{C2S / 'models' / 'amino_acid_type.ckpt'}'",
                f"atom_checkpoint: '{C2S / 'models' / 'atom_type.ckpt'}'",
                f"confidence_score_models: '{C2S / 'models'}'",
                "infer_run_on: 'gpu'",
                "infer_on_gpu: 0",
                "threshold: 0.4",
                "clustering_threshold: 2",
                "use_sequence: 'full'",
                "std_lambda: 10",
                "",
            ]
        )
    )
    script = f"""
set -e
source {CONDA}
conda activate cryo2struct
export CUDA_VISIBLE_DEVICES=0
cd {C2S}
python3 cryo2struct.py --density_map_name {emdb}
"""
    print("Cryo2Struct", pdb, flush=True)
    with log.open("w") as handle:
        proc = subprocess.run(["bash", "-c", script], stdout=handle, stderr=subprocess.STDOUT, text=True)
    if proc.returncode == 0:
        out_dir.mkdir(parents=True, exist_ok=True)
        for pdb_file in (INPUT / emdb).glob("*cryo2struct*.pdb"):
            shutil.copy2(pdb_file, out_dir / pdb_file.name)
        for pdb_file in (INPUT / emdb).glob("*.pdb"):
            if "cryo2struct" in pdb_file.name or True:
                shutil.copy2(pdb_file, out_dir / pdb_file.name)
    else:
        print("FAIL", pdb, "see", log, flush=True)
    return proc.returncode


def main() -> int:
    INPUT.mkdir(parents=True, exist_ok=True)
    failed = []
    for pdb in FAST_FIRST:
        if pdb in ("9XON", "9TBI"):
            print("skip", pdb, "known fail/too slow", flush=True)
            failed.append(pdb)
            continue
        try:
            emdb = stage_target(pdb)
            if run_one(pdb, emdb) != 0:
                failed.append(pdb)
        except Exception as exc:
            print("FAIL", pdb, exc, flush=True)
            failed.append(pdb)
    print("failed", failed)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
