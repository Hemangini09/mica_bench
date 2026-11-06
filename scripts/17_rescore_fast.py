#!/usr/bin/env python3
"""Rescore every finished model with vectorized Cα matching. Writes scores_master.csv."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import (  # noqa: E402
    AUC_DMAX,
    AUC_DMIN,
    AUC_DSTEP,
    CA_MATCH_CUTOFF,
    DATA,
    METHODS,
    MODELS,
    PANEL,
    RESULTS,
    RUN_PDB_IDS,
)
from importlib import import_module

metrics = import_module("8_evaluate_metrics")


def candidate_pairs(model, native, max_dist: float):
    m_xyz = np.asarray([xyz for _, xyz in model], dtype=np.float64)
    n_xyz = np.asarray([xyz for _, xyz in native], dtype=np.float64)
    if len(m_xyz) == 0 or len(n_xyz) == 0:
        return np.empty(0), np.empty(0, dtype=int), np.empty(0, dtype=int)
    tree = cKDTree(n_xyz)
    lists = tree.query_ball_point(m_xyz, max_dist, workers=-1)
    ii, jj, dists = [], [], []
    for i, neighbors in enumerate(lists):
        if not neighbors:
            continue
        diff = n_xyz[neighbors] - m_xyz[i]
        dist = np.sqrt(np.einsum("ij,ij->i", diff, diff))
        ii.extend([i] * len(neighbors))
        jj.extend(neighbors)
        dists.append(dist)
    if not ii:
        return np.empty(0), np.empty(0, dtype=int), np.empty(0, dtype=int)
    return np.concatenate(dists), np.asarray(ii, dtype=np.int32), np.asarray(jj, dtype=np.int32)


def greedy_at(dists, ii, jj, cutoff, model, native):
    n_model, n_native = len(model), len(native)
    if len(dists) == 0 or n_model == 0 or n_native == 0:
        return 0, 0, 0.0
    order = np.argsort(dists, kind="mergesort")
    used_m = np.zeros(n_model, dtype=bool)
    used_n = np.zeros(n_native, dtype=bool)
    matched = 0
    identical = 0
    sum_sq = 0.0
    for k in order:
        dist = float(dists[k])
        if dist > cutoff:
            break
        i = int(ii[k])
        j = int(jj[k])
        if used_m[i] or used_n[j]:
            continue
        used_m[i] = True
        used_n[j] = True
        matched += 1
        sum_sq += dist * dist
        if model[i][0] == native[j][0]:
            identical += 1
    rmsd = (sum_sq / matched) ** 0.5 if matched else 0.0
    return matched, identical, rmsd


def auc_pr_curve(curve):
    """Area under precision–recall as the distance cutoff loosens from 0.5 to 5 Å.

    A model that is already perfect at 0.5 Å is the single point (recall=1, precision=1).
    The curve is anchored at recall 0 so that point integrates to 1, not to 0.
    """
    cleaned = []
    for recall, precision in sorted(curve, key=lambda item: (item[0], -item[1])):
        if cleaned and abs(cleaned[-1][0] - recall) < 1e-9:
            if precision > cleaned[-1][1]:
                cleaned[-1] = (recall, precision)
            continue
        cleaned.append((recall, precision))
    if not cleaned:
        return 0.0
    if cleaned[0][0] > 0:
        cleaned.insert(0, (0.0, cleaned[0][1]))
    area = 0.0
    for i in range(1, len(cleaned)):
        r0, p0 = cleaned[i - 1]
        r1, p1 = cleaned[i]
        area += (r1 - r0) * (p0 + p1) / 2.0
    return round(area, 4)


def score_atoms(model, native):
    n_model, n_native = len(model), len(native)
    dists, ii, jj = candidate_pairs(model, native, AUC_DMAX)
    if len(dists):
        order = np.argsort(dists, kind="mergesort")
        dists, ii, jj = dists[order], ii[order], jj[order]

    def at(cutoff):
        matched, identical, rmsd = greedy_at(dists, ii, jj, cutoff, model, native)
        precision = matched / n_model if n_model else 0.0
        recall = matched / n_native if n_native else 0.0
        f1 = (2 * precision * recall / (precision + recall)) if precision + recall else 0.0
        coverage = n_model / n_native if n_native else 0.0
        return {
            "precision": round(precision, 4),
            "recall": round(recall, 4),
            "f1": round(f1, 4),
            "ca_rmsd": round(rmsd, 4),
            "ca_quality": round(precision * coverage, 4),
            "sequence_match": round(identical / matched, 4) if matched else 0.0,
            "n_model_ca": n_model,
            "n_native_ca": n_native,
            "n_matched_ca": matched,
            "_pr": (recall, precision),
        }

    stats = at(CA_MATCH_CUTOFF)
    curve = []
    d = AUC_DMIN
    while d <= AUC_DMAX + 1e-9:
        curve.append(at(d)["_pr"])
        d += AUC_DSTEP
    stats["auc_pr"] = auc_pr_curve(curve)
    del stats["_pr"]
    return stats


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []
    natives = {}
    for pdb in RUN_PDB_IDS:
        native = DATA / pdb / f"{pdb}.cif"
        natives[pdb] = metrics.parse_structure(native)
        print(f"native {pdb} ca {len(natives[pdb])}", flush=True)
    for pdb in RUN_PDB_IDS:
        native_atoms = natives[pdb]
        native_path = DATA / pdb / f"{pdb}.cif"
        for method in METHODS:
            folders = [MODELS / method / pdb]
            if method == "mica":
                folders.append(MODELS / "mica_struct_free" / pdb)
            model_path = None
            for folder in folders:
                model_path = metrics.find_built_model(folder)
                if model_path is not None:
                    break
            row = {
                "pdb": pdb,
                "method": method,
                "resolution": PANEL[pdb]["resolution"],
                "status": "missing",
                "model_path": "",
            }
            if model_path is None:
                rows.append(row)
                print(f"{pdb} {method} missing", flush=True)
                continue
            model_atoms = metrics.parse_structure(model_path)
            row.update(score_atoms(model_atoms, native_atoms))
            row.update(metrics.usalign_metrics(model_path, native_path))
            row["status"] = "ok"
            row["model_path"] = str(model_path)
            rows.append(row)
            print(
                f"{pdb} {method} tm {row.get('tm_score')} f1 {row.get('f1')} "
                f"rec {row.get('recall')} seq {row.get('sequence_match')} "
                f"n {row.get('n_model_ca')}/{row.get('n_native_ca')}",
                flush=True,
            )
    out = RESULTS / "scores_master.csv"
    fieldnames = sorted({k for r in rows for k in r})
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print("wrote", out, "rows", len(rows), flush=True)


if __name__ == "__main__":
    main()
