#!/usr/bin/env python3
"""Run one fast-track lane: phenix, or cryo2struct then MICA.

9RVN and 9TBI Phenix stay stopped so they can be resumed later.
Cryo2Struct runs only on the four small maps. The larger Cryo2Struct
maps are not started.
"""

from __future__ import annotations

import importlib
import sys

q = importlib.import_module("15_queue_phenix_c2s_mica")


def phenix_lane() -> int:
    failed = []
    print("Phenix fast-track", " ".join(q.PHENIX), flush=True)
    print("deferred", "9RVN", "9TBI", flush=True)
    for pdb in q.PHENIX:
        if q.phenix(pdb) != 0:
            failed.append(pdb)
    print("phenix failed", failed, flush=True)
    return 1 if failed else 0


def cryo2struct_lane() -> int:
    failed = []
    print("Cryo2Struct fast maps", " ".join(q.CRYO2STRUCT), flush=True)
    print("deferred", "9RVN", "9WP6", "9VF5", "9XDE", "9XON", flush=True)
    for pdb in q.CRYO2STRUCT:
        if q.cryo2struct(pdb) != 0:
            failed.append(f"cryo2struct:{pdb}")
    print("starting MICA", " ".join(q.MICA), flush=True)
    for pdb in q.MICA:
        if q.mica(pdb) != 0:
            failed.append(f"mica:{pdb}")
    print("c2s-mica failed", failed, flush=True)
    return 1 if failed else 0


def main() -> int:
    lane = sys.argv[1] if len(sys.argv) > 1 else ""
    if lane == "phenix":
        return phenix_lane()
    if lane == "cryo2struct":
        return cryo2struct_lane()
    raise SystemExit("usage: 16_fasttrack.py phenix|cryo2struct")


if __name__ == "__main__":
    sys.exit(main())
