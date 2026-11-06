#!/usr/bin/env python3
"""Score predicted models vs deposited CIFs. Writes results/scores_master.csv."""

from __future__ import annotations

import csv
import math
import subprocess
from pathlib import Path

from config import (
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
    USALIGN,
)


def split_cif(line: str) -> list[str]:
    parts, token, quote = [], [], None
    for char in line.strip():
        if quote:
            if char == quote:
                quote = None
            else:
                token.append(char)
        elif char in ("'", '"'):
            quote = char
        elif char.isspace():
            if token:
                parts.append("".join(token))
                token = []
        else:
            token.append(char)
    if token:
        parts.append("".join(token))
    return parts


def parse_pdb(text: str):
    atoms, seen = [], set()
    for line in text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        if line[16:17] not in (" ", "A", "."):
            continue
        key = (line[21:22], line[22:27])
        if key in seen:
            continue
        seen.add(key)
        atoms.append(
            (line[17:20].strip(), (float(line[30:38]), float(line[38:46]), float(line[46:54])))
        )
    return atoms


def parse_cif(text: str):
    lines = text.splitlines()
    header, header_end = None, 0
    for i, line in enumerate(lines):
        if line.startswith("_atom_site."):
            header = [] if header is None else header
            header.append(line[len("_atom_site.") :].strip())
            header_end = i + 1
        elif header is not None:
            break
    if not header:
        return []
    index = {name: n for n, name in enumerate(header)}
    required = ("label_atom_id", "label_comp_id", "Cartn_x", "Cartn_y", "Cartn_z")
    if any(name not in index for name in required):
        return []
    atoms, seen = [], set()
    for line in lines[header_end:]:
        if not line.startswith("ATOM"):
            if atoms:
                break
            continue
        parts = split_cif(line)
        if len(parts) <= max(index[name] for name in required):
            continue
        if "group_PDB" in index and parts[index["group_PDB"]] != "ATOM":
            continue
        if parts[index["label_atom_id"]] != "CA":
            continue
        if "label_alt_id" in index and parts[index["label_alt_id"]] not in (".", "?", "A"):
            continue
        chain = parts[index["label_asym_id"]] if "label_asym_id" in index else ""
        seq = parts[index["label_seq_id"]] if "label_seq_id" in index else ""
        key = (chain, seq)
        if key in seen:
            continue
        seen.add(key)
        atoms.append(
            (
                parts[index["label_comp_id"]],
                (
                    float(parts[index["Cartn_x"]]),
                    float(parts[index["Cartn_y"]]),
                    float(parts[index["Cartn_z"]]),
                ),
            )
        )
    return atoms


def parse_structure(path: Path):
    text = Path(path).read_text(errors="replace")
    if "_atom_site." in text[:20000] or text.lstrip().startswith("data_"):
        atoms = parse_cif(text)
        if atoms:
            return atoms
    return parse_pdb(text)


def match_pairs(model, native, cutoff: float):
    pairs = []
    for i, (_, xyz_m) in enumerate(model):
        for j, (_, xyz_n) in enumerate(native):
            dist = math.dist(xyz_m, xyz_n)
            if dist <= cutoff:
                pairs.append((dist, i, j))
    pairs.sort()
    used_m, used_n, accepted, identical = set(), set(), [], 0
    for dist, i, j in pairs:
        if i in used_m or j in used_n:
            continue
        used_m.add(i)
        used_n.add(j)
        accepted.append(dist)
        if model[i][0] == native[j][0]:
            identical += 1
    return accepted, identical


def ca_metrics(model, native, cutoff=CA_MATCH_CUTOFF):
    accepted, identical = match_pairs(model, native, cutoff)
    matched = len(accepted)
    n_model, n_native = len(model), len(native)
    precision = matched / n_model if n_model else 0.0
    recall = matched / n_native if n_native else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if precision + recall else 0.0
    ca_rmsd = math.sqrt(sum(d * d for d in accepted) / matched) if matched else 0.0
    coverage = n_model / n_native if n_native else 0.0
    quality = precision * coverage
    return {
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
        "ca_rmsd": round(ca_rmsd, 4),
        "ca_quality": round(quality, 4),
        "sequence_match": round(identical / matched, 4) if matched else 0.0,
        "n_model_ca": n_model,
        "n_native_ca": n_native,
        "n_matched_ca": matched,
    }


def auc_pr(model, native) -> float:
    """Trapezoidal AUC of precision vs recall over distance cutoffs."""
    points = []
    d = AUC_DMIN
    while d <= AUC_DMAX + 1e-9:
        stats = ca_metrics(model, native, cutoff=d)
        points.append((stats["recall"], stats["precision"]))
        d += AUC_DSTEP
    points.sort()
    area = 0.0
    for i in range(1, len(points)):
        r0, p0 = points[i - 1]
        r1, p1 = points[i]
        area += (r1 - r0) * (p0 + p1) / 2.0
    return round(area, 4)


def usalign_metrics(model: Path, native: Path):
    result = subprocess.run(
        [
            str(USALIGN),
            str(model),
            str(native),
            "-mm",
            "1",
            "-ter",
            "1",
            "-mol",
            "prot",
            "-outfmt",
            "2",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    data = None
    for row in result.stdout.splitlines():
        if row.startswith("#") or not row.strip():
            continue
        cols = row.split("\t")
        if len(cols) >= 11:
            data = cols
            break
    if data is None:
        return {
            "tm_score": None,
            "usalign_rmsd": None,
            "sequence_identity": None,
            "aligned_ca_length": None,
            "usalign_error": (result.stderr or result.stdout)[-400:],
        }
    return {
        "tm_score": float(data[3]),
        "usalign_rmsd": float(data[4]),
        "sequence_identity": float(data[7]),
        "aligned_ca_length": int(float(data[10])),
    }


def _is_real_model(path: Path) -> bool:
    if path.suffix.lower() not in {".cif", ".pdb"}:
        return False
    if path.stat().st_size <= 1000:
        return False
    if "(init)" in path.name:
        return False
    text = str(path)
    junk = (
        "TEMP_BUILD",
        "segmented_maps",
        "dots.pdb",
        "find_helices_strands",
        "helices_strands",
        "trace_and_build/",
    )
    return not any(part in text for part in junk)


def find_built_model(folder: Path):
    if not folder.exists():
        return None
    top = [p for p in folder.iterdir() if p.is_file() and _is_real_model(p)]
    if top:
        named = [p for p in top if "all_atom" in p.name]
        if not named:
            named = [p for p in top if p.name == "map_to_model.pdb" or "ca_model" in p.name or p.suffix == ".cif"]
        return max(named or top, key=lambda p: p.stat().st_size)
    deep = [p for p in folder.rglob("*") if p.is_file() and _is_real_model(p)]
    return max(deep, key=lambda p: p.stat().st_size) if deep else None


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []
    for pdb in RUN_PDB_IDS:
        native = DATA / pdb / f"{pdb}.cif"
        if not native.exists():
            print(f"missing native {pdb}")
            continue
        native_atoms = parse_structure(native)
        for method in METHODS:
            folders = [MODELS / method / pdb]
            if method == "mica":
                folders.append(MODELS / "mica_struct_free" / pdb)
            model_path = None
            for folder in folders:
                model_path = find_built_model(folder)
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
                continue
            model_atoms = parse_structure(model_path)
            row.update(ca_metrics(model_atoms, native_atoms))
            row["auc_pr"] = auc_pr(model_atoms, native_atoms)
            row.update(usalign_metrics(model_path, native))
            row["status"] = "ok"
            row["model_path"] = str(model_path)
            rows.append(row)
            print(pdb, method, row.get("tm_score"), row.get("f1"))
    out = RESULTS / "scores_master.csv"
    fieldnames = sorted({k for r in rows for k in r})
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print("wrote", out, "rows", len(rows))


if __name__ == "__main__":
    main()
