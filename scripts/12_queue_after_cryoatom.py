#!/usr/bin/env python3
"""After CryoAtom frees the GPU, run the methods that already finish.

Order: EModelX(+AF), Cryo2Struct, phenix.map_to_model.
Skips 9TBI and MICA. Only the small maps still missing these models.
"""

from __future__ import annotations

import subprocess
import sys
import time

ORDER = ("9N5N", "33AW", "9UIT", "9S3M")


def cryoatom_running() -> bool:
    active = subprocess.run(
        ["systemctl", "--user", "is-active", "--quiet", "mica-bench-cryoatom.service"]
    ).returncode == 0
    listed = subprocess.run(
        ["pgrep", "-f", "cryoatom build"], stdout=subprocess.DEVNULL
    ).returncode == 0
    return active or listed


def wait_cryoatom() -> None:
    while cryoatom_running():
        print("waiting for CryoAtom to finish", flush=True)
        time.sleep(30)


def main() -> int:
    print("queued behind CryoAtom:", " ".join(ORDER), flush=True)
    wait_cryoatom()
    print("CryoAtom done, starting EModelX", flush=True)

    from importlib import import_module

    emx = import_module("9_run_emodelx")
    phenix = import_module("7_run_phenix_map_to_model")

    failed = []
    for pdb in ORDER:
        if emx.already(pdb):
            print("emodelx already", pdb, flush=True)
            continue
        rc = emx.run_one(pdb)
        ok = emx.already(pdb)
        print("emodelx", pdb, "rc", rc, "model", ok, flush=True)
        if not ok:
            failed.append(f"emodelx:{pdb}")

    print("starting Cryo2Struct", flush=True)
    # scipy and mrcfile live in the cryo2struct env, not this interpreter.
    c2s_code = """
import importlib, sys
c2s = importlib.import_module("4_run_cryo2struct")
failed = []
order = sys.argv[1:]
for pdb in order:
    try:
        emdb = c2s.stage_target(pdb)
        rc = c2s.run_one(pdb, emdb)
        print("cryo2struct", pdb, "rc", rc, flush=True)
        if rc != 0:
            failed.append(pdb)
    except Exception as exc:
        print("cryo2struct FAIL", pdb, exc, flush=True)
        failed.append(pdb)
print("cryo2struct failed", failed, flush=True)
sys.exit(1 if failed else 0)
"""
    proc = subprocess.run(
        [
            "bash",
            "-c",
            "source /home/hemangini/miniconda3/etc/profile.d/conda.sh && "
            "conda activate cryo2struct && python3 -c \"$1\" \"${@:2}\"",
            "bash",
            c2s_code,
            *ORDER,
        ]
    )
    if proc.returncode != 0:
        failed.append("cryo2struct")

    print("starting phenix.map_to_model", flush=True)
    for pdb in ORDER:
        try:
            rc = phenix.run_one(pdb)
            print("phenix", pdb, "rc", rc, flush=True)
            if rc != 0:
                failed.append(f"phenix:{pdb}")
        except Exception as exc:
            print("phenix FAIL", pdb, exc, flush=True)
            failed.append(f"phenix:{pdb}")

    print("after-cryoatom failed", failed, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
