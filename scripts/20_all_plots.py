#!/usr/bin/env python3
"""One page: coverage, PR, TM/F1 stems, AUC, TM–F1 quadrants, and score bars."""

from __future__ import annotations

import csv
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

ROOT = Path(__file__).resolve().parents[1]
SCORES = ROOT / "results" / "scores_master.csv"
FIG = ROOT / "figures"

SERIES = (
    ("modelangelo", "ModelAngelo", "s", "#D64545"),
    ("emodelx", "EModelX(+AF)", "o", "#3A6FBF"),
    ("cryoatom", "CryoAtom", "^", "#2E8B57"),
    ("mica", "MICA", "D", "#6B4C9A"),
)
ORDER = tuple(item[0] for item in SERIES)
LABELS = {key: name for key, name, _m, _c in SERIES}
COLORS = {key: color for key, _n, _m, color in SERIES}
INK = "#1C1917"
FS = 14
FS_LETTER = 16
TICKS = (0.0, 0.25, 0.50, 0.75, 1.00)
YTICKS_COV = (0.05, 0.5, 1, 10, 60)
METRICS = (
    ("tm_score", "TM-score"),
    ("f1", "F1 at 3 Å"),
    ("precision", "Precision"),
    ("recall", "Recall"),
    ("sequence_match", "Sequence match"),
    ("auc_pr", "AUC-PR"),
)
PR_CALLOUTS = {
    ("modelangelo", "9TBI"): (0.70, 0.22),
    ("cryoatom", "9TBI"): (0.70, 0.14),
    ("cryoatom", "9RVN"): (0.70, 0.06),
    ("modelangelo", "9XON"): (0.78, 0.58),
    ("cryoatom", "9XON"): (0.78, 0.48),
    ("cryoatom", "33AW"): (0.78, 0.30),
}
PR_SHIFTS = {("modelangelo", "9RVN"): (8, 8), ("cryoatom", "9WSD"): (6, 8)}
AUC_SHIFTS = {
    ("emodelx", "9W22"): (0, -15, "center"),
    ("mica", "9W22"): (0, -29, "center"),
    ("modelangelo", "33AW"): (10, 8, "left"),
    ("cryoatom", "9WSD"): (-8, 12, "right"),
}
CORNERS = {
    ("modelangelo", "9TBI"): (0.68, 0.30),
    ("cryoatom", "9TBI"): (0.68, 0.22),
    ("modelangelo", "9RVN"): (0.68, 0.14),
    ("cryoatom", "9RVN"): (0.68, 0.06),
}
SHIFTS = {
    ("modelangelo", "9UIT"): (8, 2),
    ("modelangelo", "9S3M"): (8, -12),
    ("modelangelo", "9W22"): (-30, 6),
    ("modelangelo", "33AW"): (6, 3),
    ("modelangelo", "9WSD"): (6, 3),
    ("modelangelo", "9XON"): (6, 2),
    ("cryoatom", "9XON"): (6, -8),
    ("cryoatom", "33AW"): (-32, 8),
    ("cryoatom", "9WSD"): (6, 6),
    ("mica", "33AW"): (6, -8),
}
PITCH = 1.7
BAR = 0.48


def load():
    rows = []
    with SCORES.open() as handle:
        for row in csv.DictReader(handle):
            if row["method"] not in ORDER or row.get("status") != "ok":
                continue
            if "/backup/" in row.get("model_path", ""):
                continue
            try:
                native = float(row["n_native_ca"]) if row.get("n_native_ca") else None
                built = float(row["n_model_ca"]) if row.get("n_model_ca") else None
                parsed = {
                    "method": row["method"],
                    "pdb": row["pdb"],
                    "resolution": float(row["resolution"]) if row["resolution"] else None,
                    "tm_score": float(row["tm_score"]) if row["tm_score"] else None,
                    "f1": float(row["f1"]) if row["f1"] else None,
                    "precision": float(row["precision"]) if row["precision"] else None,
                    "recall": float(row["recall"]) if row["recall"] else None,
                    "sequence_match": float(row["sequence_match"]) if row["sequence_match"] else None,
                    "auc_pr": float(row["auc_pr"]) if row["auc_pr"] else None,
                    "coverage": (built / native) if native and native > 0 and built is not None else None,
                }
            except ValueError:
                continue
            rows.append(parsed)
    return rows


def maps_in_order(rows):
    maps, seen = [], set()
    for row in sorted(rows, key=lambda item: (item["resolution"] or 0, item["pdb"])):
        if row["pdb"] not in seen:
            seen.add(row["pdb"])
            maps.append(row["pdb"])
    resolution = {row["pdb"]: row["resolution"] for row in rows if row["resolution"] is not None}
    return maps, resolution


def hide_top_right(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color("black")
        ax.spines[side].set_linewidth(1.0)
    ax.tick_params(axis="both", which="major", direction="out", length=3.2, width=0.8, color="black", bottom=True, top=False, left=True, right=False, pad=2, labelsize=FS)


def square_frame(ax):
    ax.set_xlim(-0.06, 1.08)
    ax.set_ylim(-0.06, 1.08)
    ax.set_xticks(TICKS)
    ax.set_yticks(TICKS)
    ax.set_xticklabels(["0", "0.25", "0.50", "0.75", "1.00"])
    ax.set_yticklabels(["0", "0.25", "0.50", "0.75", "1.00"])
    ax.set_aspect("equal", adjustable="box")
    hide_top_right(ax)
    ax.spines["bottom"].set_bounds(0, 1)
    ax.spines["left"].set_bounds(0, 1)
    ax.spines["bottom"].set_zorder(5)
    ax.spines["left"].set_zorder(5)


def draw_coverage(axes, rows, maps, resolution):
    lookup = {(row["pdb"], row["method"]): row for row in rows}
    for ax, method in zip(axes, ORDER):
        color = COLORS[method]
        ax.set_yscale("log")
        ax.axhline(1, color="#D6D3D1", linewidth=0.9, zorder=1)
        for index, pdb in enumerate(maps):
            row = lookup.get((pdb, method))
            if row is None or row["coverage"] is None:
                continue
            ax.plot([index, index], [1, row["coverage"]], color=color, linewidth=1.6, solid_capstyle="round", zorder=2)
            ax.scatter([index], [row["coverage"]], s=36, color=color, zorder=3, linewidths=0)
        ax.set_xlim(-0.35, len(maps) - 0.65)
        ax.set_ylim(0.032, 110)
        ax.set_yticks(YTICKS_COV)
        ax.set_yticklabels(["0.05", "0.5", "1", "10", "60"])
        ax.set_xticks(range(len(maps)))
        ax.set_xticklabels([])
        hide_top_right(ax)
        ax.spines["bottom"].set_bounds(0, len(maps) - 1)
        ax.spines["left"].set_bounds(YTICKS_COV[0], YTICKS_COV[-1])
        ax.set_title(LABELS[method], loc="left", color=color, fontsize=FS, pad=1)
    axes[-1].set_xticklabels(maps, rotation=45, ha="right")
    axes[1].set_ylabel("Cα built / deposit")


def draw_tm_stems(axes, rows, maps):
    lookup = {(row["pdb"], row["method"]): row for row in rows}
    for ax, method in zip(axes, ORDER):
        color = COLORS[method]
        for index, pdb in enumerate(maps):
            row = lookup.get((pdb, method))
            if row is None or row["tm_score"] is None or row["f1"] is None:
                continue
            ax.plot([index, index], [row["tm_score"], row["f1"]], color=color, linewidth=1.6, solid_capstyle="round", zorder=2)
            ax.scatter([index], [row["tm_score"]], s=36, color=color, zorder=3, linewidths=0)
            ax.scatter([index], [row["f1"]], s=36, facecolor="white", edgecolor=color, linewidths=1.2, zorder=4)
        ax.set_xlim(-0.35, len(maps) - 0.65)
        ax.set_ylim(-0.06, 1.12)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_xticks(range(len(maps)))
        ax.set_xticklabels([])
        hide_top_right(ax)
        ax.spines["bottom"].set_bounds(0, len(maps) - 1)
        ax.spines["left"].set_bounds(0, 1)
        ax.set_title(LABELS[method], loc="left", color=color, fontsize=FS, pad=1)
    axes[-1].set_xticklabels(maps, rotation=45, ha="right")
    axes[1].set_ylabel("Score")
    axes[-1].legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=INK, markersize=7, label="TM-score"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=INK, markeredgewidth=1.2, markersize=7, label="Cα F1"),
        ],
        loc="lower left",
        frameon=False,
        fontsize=FS,
        handletextpad=0.3,
    )


def draw_pr(ax, rows):
    recall = np.linspace(0.001, 1, 400)
    for level in (0.5, 0.8, 0.9):
        precision = level * recall / (2 * recall - level)
        keep = (precision >= 0) & (precision <= 1) & (2 * recall > level)
        ax.plot(recall[keep], precision[keep], color="#C8C2BA", linewidth=1.15, linestyle=(0, (5, 3.5)), zorder=1)
    handles = []
    for method, label, marker, color in SERIES:
        picked = [row for row in rows if row["method"] == method]
        ax.scatter(
            [row["recall"] for row in picked],
            [row["precision"] for row in picked],
            s=52,
            marker=marker,
            color=color,
            linewidths=0.4,
            edgecolors="white",
            zorder=3,
            clip_on=False,
        )
        handles.append(
            Line2D([0], [0], marker=marker, color="none", markerfacecolor=color, markeredgecolor=color, markersize=9, label=label)
        )
    square_frame(ax)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.legend(handles=handles, loc="center left", frameon=False, fontsize=FS, handletextpad=0.35, borderaxespad=0.6)


def draw_pr_labeled(ax, rows):
    recall = np.linspace(0.001, 1, 400)
    for level, text_at in ((0.5, (0.40, 0.70)), (0.8, (0.52, 0.86)), (0.9, (0.66, 1.01))):
        precision = level * recall / (2 * recall - level)
        keep = (precision >= 0) & (precision <= 1) & (2 * recall > level)
        ax.plot(recall[keep], precision[keep], color="#C8C2BA", linewidth=1.0, zorder=1)
        ax.text(*text_at, f"F1 {level:.1f}", fontsize=FS, color="#8A847C", ha="center", va="bottom", zorder=2)
    for method, _label, marker, color in SERIES:
        picked = [row for row in rows if row["method"] == method]
        ax.scatter(
            [row["recall"] for row in picked],
            [row["precision"] for row in picked],
            s=52,
            marker=marker,
            color=color,
            linewidths=0.4,
            edgecolors="white",
            zorder=3,
            clip_on=False,
        )
        for row in picked:
            if row["precision"] is None or row["recall"] is None:
                continue
            if not (row["precision"] < 0.55 or row["recall"] < 0.45 or (row["f1"] or 1) < 0.45):
                continue
            spot = PR_CALLOUTS.get((method, row["pdb"]))
            if spot is not None:
                ax.annotate(
                    row["pdb"],
                    (row["recall"], row["precision"]),
                    xytext=spot,
                    fontsize=FS,
                    color=color,
                    ha="right",
                    va="center",
                    zorder=4,
                    arrowprops={"arrowstyle": "-", "color": color, "lw": 0.7, "shrinkA": 1, "shrinkB": 4},
                )
                continue
            ax.annotate(
                row["pdb"],
                (row["recall"], row["precision"]),
                textcoords="offset points",
                xytext=PR_SHIFTS.get((method, row["pdb"]), (5, 3)),
                fontsize=FS,
                color=color,
                zorder=4,
            )
    square_frame(ax)
    ax.set_xlabel("Recall at 3 Å, fraction of the deposit found")
    ax.set_ylabel("Precision at 3 Å, fraction of the written atoms that are real")


def draw_auc(ax, rows):
    ax.plot([0, 1], [0, 1], color="#C4C1BE", linewidth=1.05, zorder=1)
    ax.text(0.52, 0.12, "Below the line, the 3 Å match\nslips when the cutoff is tighter", fontsize=FS, color="#78716C", ha="left", va="center", zorder=2)
    for method, _label, marker, color in SERIES:
        picked = [row for row in rows if row["method"] == method and row["f1"] is not None and row["auc_pr"] is not None]
        ax.scatter(
            [row["f1"] for row in picked],
            [row["auc_pr"] for row in picked],
            s=56,
            marker=marker,
            color=color,
            linewidths=0.4,
            edgecolors="white",
            zorder=4 if method == "emodelx" else 3,
            clip_on=False,
        )
        if not picked:
            continue
        loose = max(picked, key=lambda row: row["f1"] - row["auc_pr"])
        if loose["f1"] - loose["auc_pr"] < 0.30:
            continue
        dx, dy, align = AUC_SHIFTS.get((method, loose["pdb"]), (6, 4, "left"))
        ax.annotate(loose["pdb"], (loose["f1"], loose["auc_pr"]), textcoords="offset points", xytext=(dx, dy), fontsize=FS, color=color, ha=align, va="center", zorder=5)
    square_frame(ax)
    ax.set_xlabel("F1 at 3 Å")
    ax.set_ylabel("AUC-PR, from 0.5 Å to 5 Å")


def draw_quadrant(ax, rows):
    ax.add_patch(Rectangle((0.5, 0.5), 0.5, 0.5, facecolor="#E5F6EC", edgecolor="none", zorder=0))
    ax.add_patch(Rectangle((0.5, 0.0), 0.5, 0.5, facecolor="#FDECEC", edgecolor="none", zorder=0))
    ax.add_patch(Rectangle((0.0, 0.5), 0.5, 0.5, facecolor="#E7F0FA", edgecolor="none", zorder=0))
    ax.plot([0, 1], [0, 1], color="#C4C1BE", linewidth=1.05, zorder=1)
    for method, _label, marker, color in SERIES:
        picked = [row for row in rows if row["method"] == method and row["tm_score"] is not None and row["f1"] is not None]
        ax.scatter(
            [row["tm_score"] for row in picked],
            [row["f1"] for row in picked],
            s=52,
            marker=marker,
            color=color,
            linewidths=0.4,
            edgecolors="white",
            zorder=3,
            clip_on=False,
        )
        for row in picked:
            if abs(row["tm_score"] - row["f1"]) < 0.45:
                continue
            spot = CORNERS.get((method, row["pdb"]))
            if spot is not None:
                ax.annotate(
                    row["pdb"],
                    (row["tm_score"], row["f1"]),
                    xytext=spot,
                    fontsize=FS,
                    color=color,
                    ha="right",
                    va="center",
                    zorder=4,
                    arrowprops={"arrowstyle": "-", "color": color, "lw": 0.7, "shrinkA": 1, "shrinkB": 4},
                )
                continue
            ax.annotate(
                row["pdb"],
                (row["tm_score"], row["f1"]),
                textcoords="offset points",
                xytext=SHIFTS.get((method, row["pdb"]), (6, 4)),
                fontsize=FS,
                color=color,
                zorder=4,
            )
    square_frame(ax)
    ax.set_xlabel("TM after superposition")
    ax.set_ylabel("Cα F1 at 3 Å")


def draw_bars(axes, rows):
    centers = [(len(ORDER) - 1 - index) * PITCH for index in range(len(ORDER))]
    for ax, (key, name) in zip(axes, METRICS):
        for index, method in enumerate(ORDER):
            picked = [row for row in rows if row["method"] == method and row[key] is not None]
            if not picked:
                continue
            center = centers[index]
            med = float(median(row[key] for row in picked))
            worst_row = min(picked, key=lambda row: (row[key], row["pdb"]))
            worst = float(worst_row[key])
            color = COLORS[method]
            ax.barh(center + 0.30, med, height=BAR, color=color, linewidth=0, zorder=3)
            ax.barh(center - 0.30, worst, height=BAR, color=color, alpha=0.38, linewidth=0, zorder=3)
            label_x = max(med, worst, 0) + 0.04
            ax.text(label_x, center + 0.30, f"{med:.2f}", va="center", ha="left", fontsize=FS, color=INK, clip_on=False)
            low = f"{worst:.2f}   {worst_row['pdb']}" if med - worst >= 0.15 else f"{worst:.2f}"
            ax.text(label_x, center - 0.30, low, va="center", ha="left", fontsize=FS, color=INK, clip_on=False)
        ax.set_xlim(-0.02, 1.72)
        ax.set_ylim(centers[-1] - 0.85, centers[0] + 0.85)
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(["0", "0.50", "1"])
        ax.set_yticks(centers)
        ax.set_yticklabels([LABELS[method] for method in ORDER])
        for label, method in zip(ax.get_yticklabels(), ORDER):
            label.set_color(COLORS[method])
        ax.set_title(name, loc="center", fontsize=FS, color=INK, pad=5)
        hide_top_right(ax)
        ax.spines["bottom"].set_bounds(0, 1)
        ax.spines["left"].set_bounds(centers[-1], centers[0])


def letter(ax, text, x=-0.12, y=1.08):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=FS_LETTER, fontweight="bold", color=INK, ha="left", va="bottom")


def draw(rows):
    plt.rcParams.update(
        {
            "font.family": "Liberation Sans",
            "font.size": FS,
            "axes.labelsize": FS,
            "xtick.labelsize": FS,
            "ytick.labelsize": FS,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "text.color": INK,
            "axes.labelcolor": INK,
        }
    )
    maps, resolution = maps_in_order(rows)
    fig = plt.figure(figsize=(16.6, 19.6))
    fig.suptitle(
        "Fold match and map-frame Cα accuracy of four cryo-EM model builders",
        x=0.09,
        y=0.988,
        ha="left",
        fontsize=FS + 2,
        color=INK,
    )
    outer = GridSpec(
        3,
        1,
        figure=fig,
        height_ratios=[1.05, 1.0, 1.28],
        left=0.09,
        right=0.985,
        top=0.935,
        bottom=0.035,
        hspace=0.28,
    )
    row0 = outer[0].subgridspec(1, 2, wspace=0.22)
    row1 = outer[1].subgridspec(1, 2, wspace=0.22)
    row2 = outer[2].subgridspec(1, 2, width_ratios=[1.18, 1.0], wspace=0.24)
    tm_grid = row0[0, 0].subgridspec(4, 1, hspace=0.34)
    cov_grid = row0[0, 1].subgridspec(4, 1, hspace=0.34)
    bar_grid = row2[0, 0].subgridspec(3, 2, wspace=0.42, hspace=0.55)

    tm_axes = [fig.add_subplot(tm_grid[i]) for i in range(4)]
    cov_axes = [fig.add_subplot(cov_grid[i]) for i in range(4)]
    ax_quad = fig.add_subplot(row1[0, 0])
    ax_auc = fig.add_subplot(row1[0, 1])
    bar_axes = [fig.add_subplot(bar_grid[r, c]) for r in range(3) for c in range(2)]
    ax_pr = fig.add_subplot(row2[0, 1])

    draw_tm_stems(tm_axes, rows, maps)
    draw_coverage(cov_axes, rows, maps, resolution)
    draw_quadrant(ax_quad, rows)
    draw_auc(ax_auc, rows)
    draw_bars(bar_axes, rows)
    draw_pr_labeled(ax_pr, rows)

    method_handles = [
        Line2D([0], [0], marker=marker, color="none", markerfacecolor=color, markeredgecolor=color, markersize=9, label=label)
        for _m, label, marker, color in SERIES
    ]
    fig.legend(
        handles=method_handles,
        loc="upper left",
        bbox_to_anchor=(0.52, 0.972),
        ncol=4,
        frameon=False,
        fontsize=FS,
        handletextpad=0.35,
        columnspacing=1.1,
    )
    ax_quad.legend(
        handles=[
            Patch(facecolor="#E5F6EC", edgecolor="#9FCBB0", label="Both good"),
            Patch(facecolor="#E7F0FA", edgecolor="#A9C3DE", label="Atoms in map, fold in pieces"),
            Patch(facecolor="#FDECEC", edgecolor="#E4B4B4", label="Fold right, atoms not"),
        ],
        loc="lower left",
        frameon=False,
        fontsize=FS,
        handlelength=1.2,
        handletextpad=0.35,
        borderaxespad=0.2,
    )
    bar_axes[1].legend(
        handles=[
            Patch(facecolor="#44403C", edgecolor="none", label="Dark bar, middle"),
            Patch(facecolor="#44403C", alpha=0.38, edgecolor="none", label="Light bar, worst"),
        ],
        loc="upper left",
        bbox_to_anchor=(0.0, 1.38),
        frameon=False,
        fontsize=FS,
        handlelength=1.3,
        handletextpad=0.4,
    )

    letter(tm_axes[0], "a", x=-0.08, y=1.18)
    letter(cov_axes[0], "b", x=-0.08, y=1.18)
    letter(ax_quad, "c", x=-0.14, y=1.04)
    letter(ax_auc, "d", x=-0.14, y=1.04)
    letter(bar_axes[0], "e", x=-0.22, y=1.18)
    letter(ax_pr, "f", x=-0.14, y=1.04)

    png = FIG / "all_in_one.png"
    pdf = FIG / "all_in_one.pdf"
    fig.savefig(png, dpi=180)
    fig.savefig(pdf)
    plt.close(fig)
    print(png)
    print(pdf)


if __name__ == "__main__":
    draw(load())
