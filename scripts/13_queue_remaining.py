#!/usr/bin/env python3
"""Run the methods that finish on the remaining B2 maps.

The 12 downloaded maps are the panel. 9TBI is excluded.
9N5N, 33AW, 9UIT, and 9S3M are handled by 12_queue_after_cryoatom.py.
This waits for that job, then runs CryoAtom, EModelX, Cryo2Struct, and
phenix.map_to_model on the four larger maps still missing those models.
"""

from __future__ import annotations

import subprocess
import sys
import time
from importlib import import_module

# Smallest map first. 9TBI stays out.
ORDER = ("9RVN", "9WP6", "9VF5", "9XDE")
CONDA = "source /home/hemangini/miniconda3/etc/profile.d/conda.sh"


def earlier_still_running() -> bool:
    if (
        subprocess.run(
            ["systemctl", "--user", "is-active", "--quiet", "mica-bench-after-cryoatom.service"]
        ).returncode
        == 0
    ):
        return True
    # Match the worker binaries, not this queue script.
    for pattern in (
        r"protocol=temp_flex",
        r"bin/cryoatom build",
        r"cryo2struct.py",
        r"phenix.map_to_model ",
    ):
        if subprocess.run(["pgrep", "-f", pattern], stdout=subprocess.DEVNULL).returncode == 0:
            return True
    return False


def wait_turn() -> None:
    while earlier_still_running():
        print("waiting for the four smaller maps to finish", flush=True)
        time.sleep(30)


def conda_python(env: str, code: str, extra_env: str = "") -> int:
    return subprocess.run(
        [
            "bash",
            "-c",
            f"{CONDA} && conda activate {env} && {extra_env} python3 -c \"$1\" \"${{@:2}}\"",
            "bash",
            code,
            *ORDER,
        ]
    ).returncode


def main() -> int:
    print("queued after the small maps:", " ".join(ORDER), flush=True)
    wait_turn()
    failed = []

    print("starting CryoAtom", flush=True)
    rc = conda_python(
        "CryoAtom2",
        """
import importlib, sys
mod = importlib.import_module("11_run_cryoatom_rest")
failed = []
for pdb in sys.argv[1:]:
    if mod.has_model(pdb):
        print("cryoatom already", pdb, flush=True)
        continue
    rc = mod.run_one(pdb)
    ok = mod.has_model(pdb)
    print("cryoatom", pdb, "rc", rc, "model", ok, flush=True)
    if not ok:
        failed.append(pdb)
print("cryoatom failed", failed, flush=True)
sys.exit(1 if failed else 0)
""",
        extra_env="export MPLBACKEND=Agg &&",
    )
    if rc != 0:
        failed.append("cryoatom")

    print("starting EModelX", flush=True)
    emx = import_module("9_run_emodelx")
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
    rc = conda_python(
        "cryo2struct",
        """
import importlib, sys
c2s = importlib.import_module("4_run_cryo2struct")
failed = []
for pdb in sys.argv[1:]:
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
""",
    )
    if rc != 0:
        failed.append("cryo2struct")

    print("starting phenix.map_to_model", flush=True)
    phenix = import_module("7_run_phenix_map_to_model")
    for pdb in ORDER:
        try:
            rc = phenix.run_one(pdb)
            print("phenix", pdb, "rc", rc, flush=True)
            if rc != 0:
                failed.append(f"phenix:{pdb}")
        except Exception as exc:
            print("phenix FAIL", pdb, exc, flush=True)
            failed.append(f"phenix:{pdb}")

    print("remaining failed", failed, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
