#!/usr/bin/env python3
"""Replace auc_pr with the anchored precision–recall area. Leaves every other column."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import DATA, METHODS, MODELS, RESULTS, RUN_PDB_IDS  # noqa: E402
from importlib import import_module

metrics = import_module("8_evaluate_metrics")
fast = import_module("17_rescore_fast")


def main():
    path = RESULTS / "scores_master.csv"
    rows = list(csv.DictReader(path.open()))
    natives = {}
    for pdb in RUN_PDB_IDS:
        natives[pdb] = metrics.parse_structure(DATA / pdb / f"{pdb}.cif")
        print(f"native {pdb} {len(natives[pdb])}", flush=True)
    for row in rows:
        if row["status"] not in {"ok", "partial"} or not row.get("model_path"):
            row["auc_pr"] = ""
            continue
        model_path = Path(row["model_path"])
        if not model_path.exists():
            print("gone", row["pdb"], row["method"], flush=True)
            continue
        atoms = metrics.parse_structure(model_path)
        scored = fast.score_atoms(atoms, natives[row["pdb"]])
        row["auc_pr"] = scored["auc_pr"]
        print(
            f"{row['pdb']} {row['method']} auc {row['auc_pr']} f1 {scored['f1']}",
            flush=True,
        )
    fields = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print("updated", path, flush=True)


if __name__ == "__main__":
    main()
