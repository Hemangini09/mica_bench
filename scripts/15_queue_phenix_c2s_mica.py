#!/usr/bin/env python3
"""After CryoAtom and EModelX: Phenix, then Cryo2Struct, then MICA.

MICA runs only on the small maps that look like the two it already finished
(9XON, 9W22). 9WSD wrote an empty model. 9TBI ran out of memory. The maps
above about 300 MB are left out for the same reason.

9TBI Phenix never wrote a final model. It runs after MICA so it cannot
block Cryo2Struct.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time

ROOT = "/home/hemangini/Documents/temp/mica_bench"
CONDA = "source /home/hemangini/miniconda3/etc/profile.d/conda.sh"

# 9N5N, 33AW, 9UIT, and 9S3M already have models. 9RVN and 9TBI are deferred:
# 9RVN sat in helix search for hours, and 9TBI never wrote a final model.
PHENIX = ("9WP6", "9VF5", "9XDE")
# Small maps only. 9W22, 9WSD, and 9TBI already have models.
# 9RVN, 9WP6, 9VF5, 9XDE, and 9XON are left for a later run.
CRYO2STRUCT = ("9N5N", "33AW", "9UIT", "9S3M")
# Same size class as 9XON (111 res) and 9W22 (358 res): maps under 100 MB.
MICA = ("9N5N", "33AW", "9UIT", "9S3M")


def still_busy() -> bool:
    if (
        subprocess.run(
            ["systemctl", "--user", "is-active", "--quiet", "mica-bench-finish-cryo-emx.service"]
        ).returncode
        == 0
    ):
        return True
    for pattern in (r"bin/cryoatom build", r"protocol=temp_flex"):
        if subprocess.run(["pgrep", "-f", pattern], stdout=subprocess.DEVNULL).returncode == 0:
            return True
    return False


def run_unit(name: str, inner: str, memory: str) -> int:
    unit = f"mica-bench-{name}.service"
    subprocess.run(["systemctl", "--user", "reset-failed", unit], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    log = f"{ROOT}/logs/{name}.log"
    print("start", name, flush=True)
    proc = subprocess.run(
        [
            "systemd-run",
            "--user",
            "--wait",
            "--collect",
            "--unit",
            unit,
            "--setenv",
            f"HOME={os.environ.get('HOME', '')}",
            "--property",
            f"WorkingDirectory={ROOT}",
            "--property",
            f"MemoryMax={memory}",
            "--property",
            "Restart=no",
            "/bin/bash",
            "-c",
            inner + f" >>{log} 2>&1",
        ]
    )
    print("done", name, "rc", proc.returncode, flush=True)
    return proc.returncode


def write_runner(name: str, body: str) -> str:
    folder = f"{ROOT}/logs/runners"
    os.makedirs(folder, exist_ok=True)
    path = f"{folder}/{name}.py"
    with open(path, "w") as handle:
        handle.write(body)
    return path


def phenix(pdb: str) -> int:
    path = write_runner(
        f"phenix_{pdb}",
        "import sys\n"
        f"sys.path.insert(0, {ROOT!r})\n"
        "import importlib\n"
        "ph = importlib.import_module('7_run_phenix_map_to_model')\n"
        f"rc = ph.run_one({pdb!r})\n"
        "raise SystemExit(rc)\n",
    )
    return run_unit(f"phenix-{pdb}", f"cd {ROOT} && python3 {path}", "120G")


def cryo2struct(pdb: str) -> int:
    path = write_runner(
        f"c2s_{pdb}",
        "import sys\n"
        f"sys.path.insert(0, {ROOT!r})\n"
        "import importlib\n"
        "c2s = importlib.import_module('4_run_cryo2struct')\n"
        f"pdb = {pdb!r}\n"
        "emdb = c2s.stage_target(pdb)\n"
        "rc = c2s.run_one(pdb, emdb)\n"
        "print('cryo2struct', pdb, 'rc', rc, flush=True)\n"
        "raise SystemExit(rc)\n",
    )
    inner = f"{CONDA} && conda activate cryo2struct && cd {ROOT} && python3 {path}"
    return run_unit(f"c2s-{pdb}", inner, "80G")


def mica(pdb: str) -> int:
    inner = f"cd {ROOT} && python3 6_run_mica.py --only {pdb}"
    return run_unit(f"mica-{pdb}", inner, "160G")


def main() -> int:
    print("waiting for CryoAtom and EModelX", flush=True)
    while still_busy():
        print("waiting for CryoAtom and EModelX", flush=True)
        time.sleep(30)
    failed = []
    print("starting Phenix", " ".join(PHENIX), flush=True)
    for pdb in PHENIX:
        if phenix(pdb) != 0:
            failed.append(f"phenix:{pdb}")
    print("starting Cryo2Struct", " ".join(CRYO2STRUCT), flush=True)
    for pdb in CRYO2STRUCT:
        if cryo2struct(pdb) != 0:
            failed.append(f"cryo2struct:{pdb}")
    print("starting MICA", " ".join(MICA), flush=True)
    for pdb in MICA:
        if mica(pdb) != 0:
            failed.append(f"mica:{pdb}")
    print("later failed", failed, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
