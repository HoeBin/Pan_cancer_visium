#!/usr/bin/env python3
"""Fig 3b and Extended Data Fig 7b: whole-slide median Jaccard per major-lineage pair (21
pairs) across cancer types.

Input : Source_Data_Fig3.xlsx sheet Fig.3b (cfg.read_source key Fig3b); Source_Data_Extended_Fig7.xlsx sheet
        ED Fig.7b (ExtFig7b; BRCA-downsampled cohort). One row per lineage pair x cancer type.
Output: Plots/<fig3b|extfig7b>_global_pairs.pdf, Tables/source_data_<fig3b|extfig7b>.csv under
        Output/Figure3 or Output/ExtendedFigure7.
Nothing is recomputed; --min-samples filters on the n_slides column (default 5 drops no row). Cancer order and
marker offsets come from sorted(cancer_type). Plotting code is that of
Analysis/Run_Fig3_CoEnrichment/04_Fig3b.py, so the PDFs match. --panel fig3b|extfig7b|all.
"""
from __future__ import annotations

# "Libraries and paths" ----
import argparse
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # Code/ for `common`
from common import config as cfg

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.family": "Arial", "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})

LINEAGE_ORDER = cfg.LINEAGE_ORDER
LINEAGE_COLORS = cfg.LINEAGE_COLORS
# (marker, color, open?) per cancer type; True = open marker
CANCER_STYLE = cfg.CANCER_STYLE

# "Inputs / outputs per panel" ----
# panel -> output figure key, Source Data key, file-name prefix
PANEL_SETS = {
    "fig3b":    {"figure": "Fig3",    "sheet": "Fig3b",    "prefix": "fig3b"},
    "extfig7b": {"figure": "ExtFig7", "sheet": "ExtFig7b", "prefix": "extfig7b"},
}


def draw(d: pd.DataFrame, outdir: Path, source_dir: Path, prefix: str, n_label: int = 3) -> None:
    """d: Source Data table (cell_type_1, cell_type_2, cancer_type, median_jaccard, n_slides)."""
    pairs = [(a, b) for a, b in combinations(LINEAGE_ORDER, 2)]
    cancers = sorted(d.cancer_type.unique())

    def pair_rows(a: str, b: str) -> pd.DataFrame:
        sel = ((d.cell_type_1 == a) & (d.cell_type_2 == b)) | ((d.cell_type_1 == b) & (d.cell_type_2 == a))
        return d[sel]

    d[["cell_type_1", "cell_type_2", "cancer_type", "median_jaccard", "n_slides"]].to_csv(
        source_dir / f"source_data_{prefix}.csv", index=False)
    offs = dict(zip(cancers, np.linspace(-0.28, 0.28, len(cancers))))

    fig, (ax, axp) = plt.subplots(2, 1, figsize=(12.5, 6.4), sharex=True,
                                  gridspec_kw={"height_ratios": [3.2, 1.0], "hspace": 0.06})

    for i, (a, b) in enumerate(pairs):
        g = pair_rows(a, b).set_index("cancer_type")["median_jaccard"]
        ax.boxplot([g.to_numpy()], positions=[i], widths=0.6, showfliers=False,
                   boxprops={"color": "#333", "linewidth": 1.0},
                   medianprops={"color": "#111", "linewidth": 1.4},
                   whiskerprops={"color": "#555"}, capprops={"color": "#555"}, zorder=1)
        for c, v in g.items():
            m, col, open_ = CANCER_STYLE[c]
            kw = {"marker": m, "s": 34, "linewidths": 1.2, "zorder": 3}
            if open_:
                ax.scatter(i + offs[c], v, facecolors="none", edgecolors=col, **kw)
            else:
                ax.scatter(i + offs[c], v, c=col, **kw)
        srt = g.sort_values()
        picked = list(srt.index[:n_label]) + list(srt.index[-n_label:])
        col_lab = [[i + offs[c], float(g[c]), c, g[c] >= srt.iloc[-n_label]] for c in picked]
        # within-column collision avoidance: top labels stack upward, bottom downward
        top = sorted([l for l in col_lab if l[3]], key=lambda l: l[1])
        bot = sorted([l for l in col_lab if not l[3]], key=lambda l: -l[1])
        gap = 0.068
        y_prev = None
        for x, y, c, _ in top:
            y_lab = y + 0.030 if y_prev is None else max(y + 0.030, y_prev + gap)
            ax.text(x, y_lab, c, ha="center", va="bottom", fontsize=11, color="#222", zorder=4)
            y_prev = y_lab
        y_prev = None
        for x, y, c, _ in bot:
            y_lab = y - 0.030 if y_prev is None else min(y - 0.030, y_prev - gap)
            ax.text(x, y_lab, c, ha="center", va="top", fontsize=11, color="#222", zorder=4)
            y_prev = y_lab

        # lower track: lineage-pair dumbbell
        ya, yb = 7 - LINEAGE_ORDER.index(a), 7 - LINEAGE_ORDER.index(b)
        axp.plot([i, i], [ya, yb], color="#8A8A8A", lw=4.5, solid_capstyle="round", zorder=1)
        axp.scatter([i, i], [ya, yb], c=[LINEAGE_COLORS[a], LINEAGE_COLORS[b]], s=95, zorder=2)

    ax.set_ylim(-0.22, 1.06)
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0.00", "0.25", "0.50", "0.75", "1.00"], fontsize=13)
    ax.set_ylabel("Median Jaccard score", fontsize=15)
    ax.set_xlim(-0.7, len(pairs) - 0.3)
    ax.set_xticks(range(len(pairs)))
    ax.set_xticklabels([])
    ax.tick_params(axis="x", bottom=False)
    ax.set_axisbelow(True)
    ax.grid(True, color="#DCDCDC", linewidth=0.8)

    axp.set_ylim(0.3, 7.7)
    axp.set_yticks([])
    axp.set_xticks(range(len(pairs)))
    axp.set_xticklabels([])
    axp.tick_params(axis="x", bottom=False)
    axp.set_axisbelow(True)
    axp.grid(axis="x", color="#E7E7E7", linewidth=0.7)
    axp.set_ylabel("Cell-type pair", fontsize=15)
    for y in range(1, 8):
        axp.axhline(y, color="#E7E7E7", lw=0.7, zorder=0)

    handles = []
    for c in cancers:
        m, col, open_ = CANCER_STYLE[c]
        handles.append(plt.Line2D([0], [0], marker=m, color="none", markeredgecolor=col,
                                  markerfacecolor="none" if open_ else col,
                                  markeredgewidth=1.3, markersize=7, label=c))
    ax.legend(handles=handles, title="Cancer Type", frameon=False, fontsize=11.5, title_fontsize=13.5,
              loc="center left", bbox_to_anchor=(1.005, 0.5))

    for ext in ("pdf",):
        fig.savefig(outdir / f"{prefix}_global_pairs.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("saved:", outdir / f"{prefix}_global_pairs.pdf")


def main() -> None:
    ap = argparse.ArgumentParser(description="Fig 3b / Extended Data Fig 7b panel from Source Data.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--panel", choices=["fig3b", "extfig7b", "all"], default="all",
                    help="fig3b = Source_Data_Fig3.xlsx sheet Fig.3b; extfig7b = Source_Data_Extended_Fig7.xlsx "
                         "sheet ED Fig.7b (BRCA-downsampled cohort)")
    ap.add_argument("--min-samples", type=int, default=5, help="min # slides per cancer type (n_slides column)")
    ap.add_argument("--n-label", type=int, default=3, help="label top-N and bottom-N cancers per pair")
    args = ap.parse_args()

    for name in (["fig3b", "extfig7b"] if args.panel == "all" else [args.panel]):
        spec = PANEL_SETS[name]
        outdir, source_dir = cfg.output_dirs(spec["figure"])
        d = cfg.read_source(spec["sheet"])
        d = d[d.n_slides >= args.min_samples]
        draw(d, outdir, source_dir, spec["prefix"], args.n_label)


if __name__ == "__main__":
    main()
