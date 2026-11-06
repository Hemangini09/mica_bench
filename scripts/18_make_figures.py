#!/usr/bin/env python3
"""Figures for the Cheng report. The precision–recall plane is left as its own chart."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

ROOT = Path("/home/hemangini/Documents/temp/mica_bench")
RESULTS = ROOT / "results"
LOGS = ROOT / "logs"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

COLORS = {
    "cryoatom": "#0072B2",
    "emodelx": "#009E73",
    "mica": "#D55E00",
    "modelangelo": "#CC79A7",
    "phenix_map_to_model": "#44403C",
    "cryo2struct": "#E69F00",
    "af3_alone": "#56B4E9",
}
LABELS = {
    "cryoatom": "CryoAtom",
    "emodelx": "EModelX",
    "mica": "MICA",
    "modelangelo": "ModelAngelo",
    "phenix_map_to_model": "Phenix",
    "cryo2struct": "Cryo2Struct",
    "af3_alone": "AF3 alone",
}
ORDER = [
    "cryoatom",
    "emodelx",
    "mica",
    "modelangelo",
    "phenix_map_to_model",
    "cryo2struct",
    "af3_alone",
]
MARKERS = {
    "cryoatom": "^",
    "emodelx": "D",
    "mica": "o",
    "modelangelo": "s",
    "phenix_map_to_model": "P",
    "cryo2struct": "v",
    "af3_alone": "X",
}
BINS = (
    ("1.5–2.0 Å", 0.0, 2.0),
    ("2.3–2.5 Å", 2.0, 3.0),
    ("3.1–4.0 Å", 3.0, 5.0),
)
INK = "#292524"
MUTED = "#A8A29E"
GRID = "#F5F5F4"


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.titleweight": "regular",
            "axes.labelsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 0.6,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 200,
            "pdf.fonttype": 42,
            "text.color": INK,
            "axes.labelcolor": INK,
            "xtick.color": "#57534E",
            "ytick.color": "#57534E",
        }
    )


def calm(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#E7E5E4")
    ax.spines["bottom"].set_color("#E7E5E4")
    ax.tick_params(length=0)
    ax.set_axisbelow(True)


def load():
    rows = []
    with (RESULTS / "scores_master.csv").open() as handle:
        for row in csv.DictReader(handle):
            if row.get("status") == "partial" or "/backup/" in row.get("model_path", ""):
                row["status"] = "partial"
            for key in (
                "tm_score",
                "f1",
                "precision",
                "recall",
                "sequence_match",
                "sequence_identity",
                "ca_rmsd",
                "usalign_rmsd",
                "resolution",
                "auc_pr",
                "n_model_ca",
                "n_native_ca",
                "n_matched_ca",
                "aligned_ca_length",
            ):
                raw = (row.get(key) or "").strip()
                row[key] = float(raw) if raw else None
            rows.append(row)
    return rows


def maps_in_order(rows):
    seen = {}
    for row in rows:
        seen[row["pdb"]] = row["resolution"]
    return sorted(seen, key=lambda pdb: (seen[pdb], pdb))


def lookup_of(rows):
    return {(r["pdb"], r["method"]): r for r in rows}


def figure_dumbbell(rows, key_a, key_b, filename, title, name_a, name_b):
    """One row per method. A stem joins two scores on the same map."""
    maps = maps_in_order(rows)
    lookup = lookup_of(rows)
    fig, axes = plt.subplots(len(ORDER), 1, figsize=(11.6, 8.8), sharex=True)
    for ax, method in zip(axes, ORDER):
        for j, pdb in enumerate(maps):
            row = lookup[(pdb, method)]
            if row["status"] != "ok" or row[key_a] is None or row[key_b] is None:
                continue
            a, b = row[key_a], row[key_b]
            color = COLORS[method]
            ax.plot([j, j], [a, b], color=color, linewidth=1.5, solid_capstyle="round", zorder=2)
            ax.scatter([j], [a], s=26, color=color, zorder=3, linewidths=0)
            ax.scatter(
                [j], [b], s=26, facecolor="white", edgecolor=color, linewidths=1.25, zorder=4
            )
        ax.set_ylim(-0.08, 1.12)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_title(LABELS[method], loc="left", color=COLORS[method], fontsize=10, pad=1)
        ax.axhline(0.5, color=GRID, linewidth=0.8, zorder=0)
        calm(ax)
        ax.yaxis.grid(True, color=GRID, linewidth=0.7)
    axes[-1].set_xticks(range(len(maps)))
    axes[-1].set_xticklabels(
        [f"{pdb}\n{lookup[(pdb, 'cryoatom')]['resolution']:.2f}" for pdb in maps]
    )
    axes[len(ORDER) // 2].set_ylabel("score")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=INK, markersize=6, label=name_a),
        Line2D(
            [0], [0], marker="o", color="none", markerfacecolor="white",
            markeredgecolor=INK, markeredgewidth=1.1, markersize=6, label=name_b,
        ),
    ]
    fig.legend(handles=handles, loc="upper right", frameon=False, ncol=2, bbox_to_anchor=(0.98, 1.01))
    fig.suptitle(title, x=0.07, ha="left", fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(FIG / filename, bbox_inches="tight")
    plt.close(fig)


def figure_coverage(rows):
    maps = maps_in_order(rows)
    lookup = lookup_of(rows)
    fig, axes = plt.subplots(len(ORDER), 1, figsize=(11.6, 8.8), sharex=True)
    for ax, method in zip(axes, ORDER):
        for j, pdb in enumerate(maps):
            row = lookup[(pdb, method)]
            if row["status"] != "ok" or not row["n_native_ca"] or row["n_model_ca"] is None:
                continue
            coverage = row["n_model_ca"] / row["n_native_ca"]
            color = COLORS[method]
            ax.plot([j, j], [1, coverage], color=color, linewidth=1.4, zorder=2)
            ax.scatter([j], [coverage], s=26, color=color, zorder=3, linewidths=0)
        ax.set_yscale("log")
        ax.set_ylim(0.15, 120)
        ax.axhline(1, color="#D6D3D1", linewidth=0.8, zorder=0)
        ax.set_yticks([0.5, 1, 10, 60])
        ax.set_yticklabels(["0.5", "1", "10", "60"])
        ax.set_title(LABELS[method], loc="left", color=COLORS[method], fontsize=10, pad=1)
        calm(ax)
    axes[-1].set_xticks(range(len(maps)))
    axes[-1].set_xticklabels(
        [f"{pdb}\n{lookup[(pdb, 'cryoatom')]['resolution']:.2f}" for pdb in maps]
    )
    axes[len(ORDER) // 2].set_ylabel("Cα built  /  Cα in the deposit")
    fig.suptitle("How many atoms came out, relative to the deposit", x=0.07, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(FIG / "precision_recall.png", bbox_inches="tight")
    plt.close(fig)


def figure_scorecard(rows):
    metrics = [
        ("tm_score", "TM-score"),
        ("f1", "F1 at 3 Å"),
        ("precision", "Precision"),
        ("recall", "Recall"),
        ("sequence_match", "Sequence match"),
        ("auc_pr", "AUC-PR"),
    ]
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    n_metrics = len(metrics)
    for j, (key, _) in enumerate(metrics):
        for i, method in enumerate(ORDER):
            data = [
                r[key]
                for r in rows
                if r["method"] == method and r["status"] == "ok" and r[key] is not None
            ]
            if not data:
                continue
            y = (n_metrics - 1 - j) + (i - (len(ORDER) - 1) / 2) * 0.11
            med = float(np.median(data))
            worst = float(np.min(data))
            color = COLORS[method]
            ax.plot([worst, med], [y, y], color=color, linewidth=1.15, zorder=2)
            ax.scatter([worst], [y], s=16, facecolor="white", edgecolor=color, linewidths=1.05, zorder=3)
            ax.scatter([med], [y], s=32, color=color, zorder=4, linewidths=0)
    ax.set_yticks(range(n_metrics))
    ax.set_yticklabels([name for _, name in reversed(metrics)])
    ax.set_xlim(-0.03, 1.05)
    ax.set_xlabel("score")
    ax.set_title("Filled dot, median. Open dot, the worst finished map.", loc="left")
    calm(ax)
    ax.xaxis.grid(True, color=GRID, linewidth=0.7)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS[m], markersize=6, label=LABELS[m])
        for m in ORDER
    ]
    ax.legend(handles=handles, frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "scorecard.png", bbox_inches="tight")
    plt.close(fig)


def figure_agreement(rows):
    fig, ax = plt.subplots(figsize=(7.6, 6.5))
    ax.add_patch(Rectangle((0.5, 0.5), 0.58, 0.58, facecolor="#F3FAF6", edgecolor="none", zorder=0))
    ax.add_patch(Rectangle((0.5, -0.06), 0.58, 0.56, facecolor="#FDF6F4", edgecolor="none", zorder=0))
    ax.add_patch(Rectangle((-0.06, 0.5), 0.56, 0.58, facecolor="#F4F8FC", edgecolor="none", zorder=0))
    ax.plot([0, 1], [0, 1], color="#D6D3D1", linewidth=0.9, zorder=1)
    for method in ORDER:
        xs, ys = [], []
        for r in rows:
            if r["method"] != method or r["status"] != "ok":
                continue
            if r["tm_score"] is None or r["f1"] is None:
                continue
            xs.append(r["tm_score"])
            ys.append(r["f1"])
            if method != "af3_alone" and abs(r["tm_score"] - r["f1"]) >= 0.45:
                ax.annotate(
                    r["pdb"], (r["tm_score"], r["f1"]),
                    textcoords="offset points", xytext=(5, 4), fontsize=7.5, color="#44403C",
                )
        ax.scatter(
            xs, ys, s=42, marker=MARKERS[method], color=COLORS[method],
            linewidths=0.4, edgecolors="white", zorder=3, label=LABELS[method],
        )
    ax.text(0.98, 0.96, "both high", ha="right", va="top", fontsize=8, color="#3F6212")
    ax.text(0.98, 0.04, "fold only", ha="right", va="bottom", fontsize=8, color="#9F1239")
    ax.text(0.04, 0.96, "atoms only", ha="left", va="top", fontsize=8, color="#1E3A8A")
    ax.set_xlim(-0.04, 1.06)
    ax.set_ylim(-0.04, 1.06)
    ax.set_xlabel("TM-score, after superposition")
    ax.set_ylabel("Cα F1 at 3 Å, in the map")
    ax.set_title("The two scores, one point per model", loc="left")
    ax.set_aspect("equal", adjustable="box")
    calm(ax)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "tm_vs_f1.png", bbox_inches="tight")
    plt.close(fig)


def figure_auc(rows):
    fig, ax = plt.subplots(figsize=(7.6, 6.5))
    ax.plot([0, 1], [0, 1], color="#D6D3D1", linewidth=0.9, zorder=1)
    labeled = []
    for method in ORDER:
        xs, ys = [], []
        best = None
        for r in rows:
            if r["method"] != method or r["status"] != "ok":
                continue
            if r["f1"] is None or r["auc_pr"] is None:
                continue
            xs.append(r["f1"])
            ys.append(r["auc_pr"])
            gap = r["f1"] - r["auc_pr"]
            if method != "af3_alone" and gap >= 0.35 and (best is None or gap > best[0]):
                best = (gap, r["pdb"], r["f1"], r["auc_pr"])
        ax.scatter(
            xs, ys, s=42, marker=MARKERS[method], color=COLORS[method],
            linewidths=0.4, edgecolors="white", zorder=3, label=LABELS[method],
        )
        if best:
            labeled.append((best, COLORS[method]))
    for (gap, pdb, x, y), color in labeled:
        ax.annotate(
            pdb, (x, y), textcoords="offset points", xytext=(6, -9),
            fontsize=7.5, color="#44403C",
        )
    ax.set_xlim(-0.04, 1.06)
    ax.set_ylim(-0.04, 1.06)
    ax.set_xlabel("Cα F1 at 3 Å")
    ax.set_ylabel("AUC-PR, cutoffs 0.5–5 Å")
    ax.set_title("A high F1 can still be a loose model", loc="left")
    ax.set_aspect("equal", adjustable="box")
    calm(ax)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "auc_f1.png", bbox_inches="tight")
    plt.close(fig)


def figure_continuity(rows):
    fig, ax = plt.subplots(figsize=(7.6, 6.5))
    for method in ORDER:
        xs, ys = [], []
        for r in rows:
            if r["method"] != method or r["status"] != "ok":
                continue
            if not r["n_native_ca"] or r["aligned_ca_length"] is None or r["f1"] is None:
                continue
            frac = r["aligned_ca_length"] / r["n_native_ca"]
            xs.append(frac)
            ys.append(r["f1"])
            if method != "af3_alone" and frac < 0.5 and r["f1"] >= 0.6:
                ax.annotate(
                    r["pdb"], (frac, r["f1"]),
                    textcoords="offset points", xytext=(5, 4), fontsize=7.5, color="#44403C",
                )
        ax.scatter(
            xs, ys, s=42, marker=MARKERS[method], color=COLORS[method],
            linewidths=0.4, edgecolors="white", zorder=3, label=LABELS[method],
        )
    ax.axvline(1, color="#E7E5E4", linewidth=0.8, zorder=0)
    ax.set_xlim(-0.03, 1.35)
    ax.set_ylim(-0.04, 1.06)
    ax.set_xlabel("Fraction of the deposit in one US-align superposition")
    ax.set_ylabel("Cα F1 at 3 Å")
    ax.set_title("Locally right, and still one chain", loc="left")
    calm(ax)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "continuity.png", bbox_inches="tight")
    plt.close(fig)


def figure_resolution(rows):
    fig, axes = plt.subplots(2, 4, figsize=(11.8, 6.0), sharex=True, sharey=True)
    flat = list(axes.ravel())
    for ax, method in zip(flat, ORDER):
        for r in rows:
            if r["method"] != method or r["status"] != "ok":
                continue
            color = COLORS[method]
            if r["tm_score"] is not None:
                ax.scatter(r["resolution"], r["tm_score"], s=28, color=color, zorder=3, linewidths=0)
            if r["f1"] is not None:
                ax.scatter(
                    r["resolution"], r["f1"], s=28, facecolor="white",
                    edgecolor=color, linewidths=1.15, zorder=4,
                )
        ax.axhline(0.5, color="#E7E5E4", linewidth=0.8, zorder=0)
        ax.set_title(LABELS[method], loc="left", color=COLORS[method], fontsize=10)
        ax.set_xlim(1.35, 4.2)
        ax.set_ylim(-0.06, 1.1)
        calm(ax)
        ax.yaxis.grid(True, color=GRID, linewidth=0.7)
    flat[-1].axis("off")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=INK, markersize=6, label="TM-score"),
        Line2D(
            [0], [0], marker="o", color="none", markerfacecolor="white",
            markeredgecolor=INK, markeredgewidth=1.1, markersize=6, label="Cα F1",
        ),
    ]
    flat[-1].legend(handles=handles, loc="center", frameon=False, fontsize=9)
    for ax in axes[1, :3]:
        ax.set_xlabel("Resolution (Å)")
    axes[0, 0].set_ylabel("score")
    axes[1, 0].set_ylabel("score")
    fig.suptitle("Each map, without connecting different proteins", x=0.02, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG / "tm_vs_resolution.png", bbox_inches="tight")
    plt.close(fig)


def figure_prior(rows):
    lookup = lookup_of(rows)
    maps = sorted(maps_in_order(rows), key=lambda pdb: lookup[(pdb, "af3_alone")]["tm_score"] or 0)
    builders = [m for m in ORDER if m != "af3_alone"]
    fig, axes = plt.subplots(2, 3, figsize=(11.4, 7.6), sharex=True)
    ys = np.arange(len(maps))
    for ax, method in zip(axes.ravel(), builders):
        for y, pdb in enumerate(maps):
            prior = lookup[(pdb, "af3_alone")]["tm_score"]
            built = lookup[(pdb, method)]
            if prior is None:
                continue
            color = COLORS[method]
            if built["status"] == "ok" and built["tm_score"] is not None:
                ax.plot([prior, built["tm_score"]], [y, y], color="#E7E5E4", linewidth=1.3, zorder=1)
                ax.scatter([built["tm_score"]], [y], s=26, color=color, zorder=3, linewidths=0)
                if built["f1"] is not None and built["f1"] < 0.4:
                    ax.annotate(
                        f"{built['f1']:.2f}", (built["tm_score"], y),
                        textcoords="offset points", xytext=(0, 7), ha="center", fontsize=6.5, color=MUTED,
                    )
            ax.scatter([prior], [y], s=18, color="#D6D3D1", zorder=2, linewidths=0)
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{pdb}  {lookup[(pdb, 'af3_alone')]['resolution']:.2f}" for pdb in maps], fontsize=7)
        ax.set_xlim(-0.04, 1.15)
        ax.set_title(LABELS[method], loc="left", color=COLORS[method], fontsize=10)
        calm(ax)
        ax.axvline(0.5, color="#F5F5F4", linewidth=0.8, zorder=0)
    for ax in axes[:, 1:].ravel():
        ax.tick_params(labelleft=False)
    for ax in axes[:, 0]:
        ax.set_ylabel("")
    for ax in axes[1, :]:
        ax.set_xlabel("TM-score")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#D6D3D1", markersize=5, label="AlphaFold 3"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=INK, markersize=6, label="builder"),
    ]
    fig.legend(handles=handles, loc="upper right", frameon=False, ncol=2, bbox_to_anchor=(0.98, 1.02))
    fig.suptitle("Gray is the prior. Color is the builder. A number is F1 below 0.40.", x=0.02, ha="left", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG / "tm_vs_af3.png", bbox_inches="tight")
    plt.close(fig)


def parse_duration(text: str) -> float | None:
    text = text.strip()
    match = re.search(r"\((\d+(?:\.\d+)?) seconds total\)", text)
    if match:
        return float(match.group(1)) / 60.0
    match = re.search(r"Run time (\d+(?:\.\d+)?) seconds", text)
    if match:
        return float(match.group(1)) / 60.0
    match = re.search(r"Total processing time:\s*(\d+(?:\.\d+)?) seconds", text)
    if match:
        return float(match.group(1)) / 60.0
    hours = re.search(r"(\d+)h", text)
    minutes = re.search(r"(\d+)min", text)
    seconds = re.search(r"(\d+(?:\.\d+)?)s", text)
    if hours or minutes or seconds:
        return (
            (int(hours.group(1)) if hours else 0) * 60
            + (int(minutes.group(1)) if minutes else 0)
            + (float(seconds.group(1)) if seconds else 0) / 60
        )
    return None


def collect_runtimes():
    found = []
    text = (LOGS / "finish_cryo_emx.log").read_text(errors="replace")
    for block in re.split(r"(?=^start )", text, flags=re.M):
        name = re.match(r"start (\S+)", block)
        if not name or "status=0/SUCCESS" not in block:
            continue
        runtime = re.search(r"Service runtime: (.+)", block)
        if not runtime:
            continue
        kind, pdb = name.group(1).split("-", 1)
        method = {"cryo": "cryoatom", "emx": "emodelx"}.get(kind)
        minutes = parse_duration(runtime.group(1))
        if method and minutes and minutes > 0.2:
            found.append((method, pdb, minutes))
    text = (LOGS / "c2s_fast.log").read_text(errors="replace")
    for block in re.split(r"(?=^start )", text, flags=re.M):
        name = re.match(r"start (\S+)", block)
        if not name or "status=0/SUCCESS" not in block:
            continue
        runtime = re.search(r"Service runtime: (.+)", block)
        if not runtime:
            continue
        kind, pdb = name.group(1).split("-", 1)
        method = {"c2s": "cryo2struct", "mica": "mica"}.get(kind)
        minutes = parse_duration(runtime.group(1))
        if method and minutes:
            found.append((method, pdb, minutes))
    for path in LOGS.glob("map_to_model_*.log"):
        pdb = path.stem.replace("map_to_model_", "")
        if len(pdb) not in (4, 5):
            continue
        for line in path.read_text(errors="replace").splitlines():
            if "wall clock time:" in line:
                minutes = parse_duration(line)
                if minutes:
                    found.append(("phenix_map_to_model", pdb, minutes))
    for path in LOGS.glob("cryo2struct_*.log"):
        pdb = path.stem.replace("cryo2struct_", "")
        body = path.read_text(errors="replace")
        for line in body.splitlines():
            if "Cryo2Struct Alignment: Run time" in line:
                minutes = parse_duration(line)
                if minutes:
                    found.append(("cryo2struct", pdb, minutes))
    for path in LOGS.glob("mica_fast_*.log"):
        pdb = path.stem.replace("mica_fast_", "")
        body = path.read_text(errors="replace")
        for line in body.splitlines():
            if "Total processing time:" in line:
                minutes = parse_duration(line)
                if minutes:
                    found.append(("mica", pdb, minutes))
    best = {}
    for method, pdb, minutes in found:
        key = (method, pdb)
        if key not in best or minutes > best[key]:
            best[key] = minutes
    return [(method, pdb, minutes) for (method, pdb), minutes in best.items()]


def figure_runtime():
    points = collect_runtimes()
    methods = [m for m in ORDER if m != "af3_alone" and any(p[0] == m for p in points)]
    fig, ax = plt.subplots(figsize=(9.4, 4.2))
    for i, method in enumerate(methods):
        y = len(methods) - 1 - i
        vals = [minutes for m, _, minutes in points if m == method]
        color = COLORS[method]
        ax.scatter(vals, [y] * len(vals), s=36, color=color, zorder=3, linewidths=0)
        med = float(np.median(vals))
        ax.plot([med, med], [y - 0.22, y + 0.22], color=color, linewidth=1.6, zorder=2)
        slow_pdb, slow = max(((pdb, minutes) for m, pdb, minutes in points if m == method), key=lambda item: item[1])
        label = f"{slow_pdb}  {slow:.0f} min" if slow >= 10 else f"{slow_pdb}  {slow:.1f} min"
        ax.annotate(label, (slow, y), textcoords="offset points", xytext=(7, 0), va="center", fontsize=7.5, color="#44403C")
    ax.set_yticks(range(len(methods)))
    ax.set_yticklabels([LABELS[m] for m in reversed(methods)])
    ax.set_xscale("log")
    ax.set_xlabel("Wall time (minutes)")
    ax.set_title("Each dot is one finished map. The tick is the median.", loc="left")
    calm(ax)
    ax.xaxis.grid(True, which="major", color=GRID, linewidth=0.7)
    fig.tight_layout()
    fig.savefig(FIG / "runtime.png", bbox_inches="tight")
    plt.close(fig)
    return points


def figure_pr_plane(rows):
    fig, ax = plt.subplots(figsize=(7.4, 6.6))
    recalls = np.linspace(0.001, 1, 400)
    for level, text in ((0.5, "F1 0.5"), (0.8, "F1 0.8"), (0.9, "F1 0.9")):
        prec = level * recalls / (2 * recalls - level)
        ok = (prec >= 0) & (prec <= 1) & (2 * recalls > level)
        ax.plot(recalls[ok], prec[ok], color="#D9D3CB", linewidth=0.8, zorder=1)
        xs, ys = recalls[ok], prec[ok]
        if len(xs):
            i = len(xs) // 3
            ax.text(xs[i], ys[i], text, fontsize=7, color="#8C8680", ha="left", va="bottom")
    for method in ORDER:
        xs, ys = [], []
        for r in rows:
            if r["method"] != method or r["status"] != "ok" or r["recall"] is None:
                continue
            xs.append(r["recall"])
            ys.append(r["precision"])
            if r["precision"] < 0.55 or r["recall"] < 0.45 or r["f1"] < 0.45:
                if method != "af3_alone":
                    ax.annotate(
                        r["pdb"],
                        (r["recall"], r["precision"]),
                        textcoords="offset points",
                        xytext=(3, 3),
                        fontsize=7,
                        color="#3A322C",
                    )
        ax.scatter(
            xs, ys, s=36, marker=MARKERS[method], color=COLORS[method],
            linewidths=0, zorder=3, label=LABELS[method],
        )
    ax.set_xlim(-0.03, 1.05)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel("Cα recall at 3 Å   (fraction of the deposit found)")
    ax.set_ylabel("Cα precision at 3 Å   (fraction of the model that is real)")
    ax.set_title("Found the chain, or buried it in extra atoms", loc="left")
    ax.set_aspect("equal", adjustable="box")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="center right", bbox_to_anchor=(1.38, 0.55))
    fig.tight_layout()
    fig.savefig(FIG / "pr_plane.png", bbox_inches="tight")
    plt.close(fig)


def write_raw_csv(rows):
    columns = [
        "pdb", "resolution_A", "method", "status", "tm_score", "f1", "precision", "recall",
        "auc_pr", "sequence_match", "sequence_identity", "ca_rmsd_A", "usalign_rmsd_A",
        "aligned_ca", "n_model_ca", "n_native_ca", "n_matched_ca", "coverage", "model_path",
    ]
    out = RESULTS / "scores_raw.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        ordered = sorted(
            rows,
            key=lambda r: (r["resolution"] or 0, r["pdb"], ORDER.index(r["method"]) if r["method"] in ORDER else 99),
        )
        for r in ordered:
            n_model, n_native = r["n_model_ca"], r["n_native_ca"]
            coverage = round(n_model / n_native, 4) if n_model and n_native else ""

            def num(key):
                val = r.get(key)
                if val is None or val == "":
                    return ""
                if key in {"aligned_ca_length", "n_model_ca", "n_native_ca", "n_matched_ca"}:
                    return int(val)
                return round(float(val), 4)

            writer.writerow(
                {
                    "pdb": r["pdb"],
                    "resolution_A": r["resolution"],
                    "method": r["method"],
                    "status": r["status"],
                    "tm_score": num("tm_score"),
                    "f1": num("f1"),
                    "precision": num("precision"),
                    "recall": num("recall"),
                    "auc_pr": num("auc_pr"),
                    "sequence_match": num("sequence_match"),
                    "sequence_identity": num("sequence_identity"),
                    "ca_rmsd_A": num("ca_rmsd"),
                    "usalign_rmsd_A": num("usalign_rmsd"),
                    "aligned_ca": num("aligned_ca_length"),
                    "n_model_ca": num("n_model_ca"),
                    "n_native_ca": num("n_native_ca"),
                    "n_matched_ca": num("n_matched_ca"),
                    "coverage": coverage,
                    "model_path": r.get("model_path") or "",
                }
            )
    print("wrote", out)


def digest(rows):
    lines = []
    for method in ORDER:
        done = [r for r in rows if r["method"] == method and r["status"] == "ok"]
        lines.append(f"{method} n={len(done)}")
    text = "\n".join(lines) + "\n"
    (RESULTS / "metric_digest.txt").write_text(text)


def main():
    style()
    rows = load()
    figure_dumbbell(
        rows, "tm_score", "f1", "tm_heatmap.png",
        "Filled dot, TM-score. Open dot, Cα F1. The stem is the disagreement.",
        "TM-score", "Cα F1",
    )
    figure_dumbbell(
        rows, "recall", "sequence_match", "recall_sequence.png",
        "Filled dot, recall. Open dot, sequence match of those atoms.",
        "Recall", "Sequence match",
    )
    figure_coverage(rows)
    figure_scorecard(rows)
    figure_agreement(rows)
    figure_auc(rows)
    figure_continuity(rows)
    figure_resolution(rows)
    figure_prior(rows)
    figure_runtime()
    figure_pr_plane(rows)
    write_raw_csv(rows)
    digest(rows)
    print("figures written")


if __name__ == "__main__":
    main()
