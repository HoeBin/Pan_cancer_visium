#!/usr/bin/env python3
"""Extended Data Fig 5b and 7d: per-slide median Jaccard by subtype-pair category (Epi-Epi, Epi-TME, TME-TME)
compared across compartments.

Input : Source_Data_Extended_Fig5.xlsx sheet ED Fig.5b (cfg.read_source key ExtFig5b; 268-slide
        cohort); Source_Data_Extended_Fig7.xlsx sheet ED Fig.7d (ExtFig7d; BRCA-downsampled cohort).
        One row per (category, slide).
Output: Plots/<panel>_pair_categories.pdf, Tables/<panel>_statistics.csv, Tables/source_data_<panel>.csv under
        Output/ExtendedFigure5 (extfig5b) or Output/ExtendedFigure7 (extfig7d).
Only the paired two-sided Wilcoxon signed-rank tests (slides with all three compartments; p = 1 if all
differences are 0) are recomputed. Plotting and test code is that of
Analysis/Run_Fig3_CoEnrichment/05_ExtDataFig5b_pair_categories.py.
"""
from __future__ import annotations

# "Libraries and paths" ----
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # Code/ for `common`
from common import config as cfg

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

plt.rcParams.update({
    "font.family": "Arial", "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})

COMP_COLORS = cfg.COMP_COLORS
ORDER = cfg.COMP_ORDER
LABELS = cfg.COMP_LABELS
CATS = cfg.PAIR_CATEGORIES

# "Inputs / outputs per panel" ----
# panel -> output figure key, Source Data key; the panel name is also the file-name prefix
PANEL_SETS = {
    "extfig5b": {"figure": "ExtFig5", "sheet": "ExtFig5b"},
    "extfig7d": {"figure": "ExtFig7", "sheet": "ExtFig7d"},
}


def fmtp(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 1e-3 else f"{p:.3f}"


def compartment_matrix(df: pd.DataFrame, index) -> pd.DataFrame:
    """Sheet columns Malignant/Boundary/Normal -> slide x compartment table with code columns ORDER (Mal/Bdy/Normal),
    as the Analysis script builds it from the merged table."""
    m = df.set_index(index)[[LABELS[c] for c in ORDER]].copy()
    m.columns = ORDER
    return m


def draw(d: pd.DataFrame, panel: str, plots: Path, tables: Path) -> None:
    """d: Source Data table (category, sample_id, Malignant, Boundary, Normal)."""
    stats_rows, src_frames = [], []
    fig, axes = plt.subplots(1, 3, figsize=(5.8, 4.4))
    for ax, cat in zip(axes, CATS):
        m = compartment_matrix(d[d.category == cat], "sample_id").dropna(subset=ORDER)
        src = m[ORDER].rename(columns=LABELS).reset_index()
        src.insert(0, "category", cat)
        src_frames.append(src)
        rng = np.random.default_rng(0)
        for i, c in enumerate(ORDER):
            v = m[c].to_numpy()
            ax.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v, s=6, c="#999",
                       alpha=0.45, linewidths=0, zorder=1)
            ax.boxplot([v], positions=[i], widths=0.5, showfliers=False, patch_artist=True,
                       zorder=2,
                       boxprops={"facecolor": "none", "edgecolor": COMP_COLORS[c], "linewidth": 1.6},
                       medianprops={"color": COMP_COLORS[c], "linewidth": 1.8},
                       whiskerprops={"color": COMP_COLORS[c]}, capprops={"color": COMP_COLORS[c]})
        ymax = float(np.nanquantile(m.to_numpy(), 0.999))
        for (a, b), h in [((0, 1), 1.03), ((1, 2), 1.13), ((0, 2), 1.23)]:
            pair = m[[ORDER[a], ORDER[b]]].dropna()
            diff = pair.iloc[:, 0] - pair.iloc[:, 1]
            p = 1.0 if np.all(diff == 0) else float(wilcoxon(pair.iloc[:, 0], pair.iloc[:, 1]).pvalue)
            stats_rows.append({"panel": panel, "category": cat,
                               "comparison": f"{ORDER[a]}-{ORDER[b]}",
                               "test": "wilcoxon_signed_rank", "n_slides": len(pair), "pvalue": p})
            y = ymax * h
            ax.plot([a, a, b, b], [y, y * 1.01, y * 1.01, y], c="#333", lw=0.8)
            ax.text(1.0, y * 1.02, f"$p$={fmtp(p)}", ha="center", fontsize=11)
        ax.set_xticks(range(3))
        ax.set_xticklabels([LABELS[c] for c in ORDER], fontsize=13, rotation=45, ha="right")
        ax.tick_params(axis="y", labelsize=13)
        ax.set_title(cat, fontsize=16)
        if cat == CATS[0]:
            ax.set_ylabel("Median Jaccard score", fontsize=15)
    fig.tight_layout()
    fig.savefig(plots / f"{panel}_pair_categories.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(stats_rows).to_csv(tables / f"{panel}_statistics.csv", index=False)
    pd.concat(src_frames).to_csv(tables / f"source_data_{panel}.csv", index=False)
    print("saved:", plots / f"{panel}_pair_categories.pdf")


def main() -> None:
    ap = argparse.ArgumentParser(description="Extended Data Fig 5b / 7d (pair categories) from Source Data.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--panel", choices=["extfig5b", "extfig7d", "all"], default="all",
                    help="extfig5b = full cohort (Source_Data_Extended_Fig5.xlsx, ED Fig.5b); "
                         "extfig7d = BRCA-downsampled cohort (Source_Data_Extended_Fig7.xlsx, ED Fig.7d)")
    args = ap.parse_args()

    for panel in (["extfig5b", "extfig7d"] if args.panel == "all" else [args.panel]):
        spec = PANEL_SETS[panel]
        plots, tables = cfg.output_dirs(spec["figure"])
        draw(cfg.read_source(spec["sheet"]), panel, plots, tables)


if __name__ == "__main__":
    main()
