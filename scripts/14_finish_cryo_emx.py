#!/usr/bin/env python3
"""Finish CryoAtom and EModelX on every downloaded map still missing a model.

Each map is its own systemd unit, so one out-of-memory kill does not stop the rest.
9WSD CryoAtom already exists. 9WSD EModelX is retried because the last file was empty.
"""

from __future__ import annotations

import subprocess
import sys
import time

ROOT = "/home/hemangini/Documents/temp/mica_bench"
CONDA = "source /home/hemangini/miniconda3/etc/profile.d/conda.sh"
# Smallest unfinished map first. 9TBI is the 864 MB map and goes last.
MAPS = ("9RVN", "9WP6", "9VF5", "9XDE", "9TBI")


def run_unit(name: str, inner: str) -> int:
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
            f"HOME={__import__('os').environ.get('HOME', '')}",
            "--property",
            f"WorkingDirectory={ROOT}",
            "--property",
            "MemoryMax=160G",
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
    path = f"{ROOT}/logs/runners/{name}.py"
    subprocess.run(["mkdir", "-p", f"{ROOT}/logs/runners"], check=True)
    with open(path, "w") as handle:
        handle.write(body)
    return path


def cryoatom(pdb: str) -> int:
    path = write_runner(
        f"cryo_{pdb}",
        "import sys\n"
        f"sys.path.insert(0, {ROOT!r})\n"
        "import importlib\n"
        "mod = importlib.import_module('11_run_cryoatom_rest')\n"
        f"pdb = {pdb!r}\n"
        "if mod.has_model(pdb):\n"
        "    print('cryoatom already', pdb, flush=True)\n"
        "    raise SystemExit(0)\n"
        "rc = mod.run_one(pdb)\n"
        "ok = mod.has_model(pdb)\n"
        "print('cryoatom', pdb, 'rc', rc, 'model', ok, flush=True)\n"
        "raise SystemExit(0 if ok else 1)\n",
    )
    inner = (
        f"{CONDA} && conda activate CryoAtom2 && export MPLBACKEND=Agg && "
        f"cd {ROOT} && python3 {path}"
    )
    return run_unit(f"cryo-{pdb}", inner)


def emodelx(pdb: str) -> int:
    path = write_runner(
        f"emx_{pdb}",
        "import sys\n"
        f"sys.path.insert(0, {ROOT!r})\n"
        "import importlib\n"
        "emx = importlib.import_module('9_run_emodelx')\n"
        f"pdb = {pdb!r}\n"
        "if emx.already(pdb):\n"
        "    print('emodelx already', pdb, flush=True)\n"
        "    raise SystemExit(0)\n"
        "rc = emx.run_one(pdb)\n"
        "ok = emx.already(pdb)\n"
        "print('emodelx', pdb, 'rc', rc, 'model', ok, flush=True)\n"
        "raise SystemExit(0 if ok else 1)\n",
    )
    inner = f"cd {ROOT} && python3 {path}"
    return run_unit(f"emx-{pdb}", inner)


def main() -> int:
    failed = []
    print("CryoAtom then EModelX:", " ".join(MAPS), "+ EModelX 9WSD", flush=True)
    for pdb in MAPS:
        if cryoatom(pdb) != 0:
            failed.append(f"cryo:{pdb}")
        # Let the GPU release before the next method.
        time.sleep(5)
        if emodelx(pdb) != 0:
            failed.append(f"emx:{pdb}")
        time.sleep(5)
    if emodelx("9WSD") != 0:
        failed.append("emx:9WSD")
    print("finish failed", failed, flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
