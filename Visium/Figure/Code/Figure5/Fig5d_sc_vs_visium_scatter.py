#!/usr/bin/env python3
"""Fig. 5d: per-cancer scatter of Endothelial-Mural LR pairs significant in both
scRNA-seq (CellPhoneDB) and Visium.

Input : Source_Data_Fig5.xlsx, cfg.read_source "Fig5d_left" (every plotted point) and
        "Fig5d_right" (per-cancer statistics).
Output: Plots/sc_vs_visium_bothsig_scatter_by_cancer_q95scaled_platformlists.pdf;
        Tables/source_data_scatter_points.csv and source_data_per_cancer_spearman.csv (both sheets as
        read) under Output/Figure5/.
Axes are x_sc_scaled / y_visium_scaled as in the sheet (platform score over the 95th percentile of its
full significant list, capped at 1), with KDE bands and y = x. The forest plot takes rho, 95% CI and n
from the right sheet as is; rho is recomputed only as a printed self-check. The 4 x 4 grid needs exactly
13 cancer types. Plotting code is that of
Analysis/Run_Fig5d_Statistics/03_Fig5d_scatter_q95_platformlists.py, so the PDF matches.
"""
from __future__ import annotations

# "Libraries and paths" ----
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # Code/ for `common`
from common import config as cfg

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm
import matplotlib.patheffects as pe
from scipy.stats import spearmanr

# "Inputs" ----
# Source_Data_Fig5.xlsx 'Fig.5d left' / 'Fig.5d right'
LEFT_KEY, RIGHT_KEY = "Fig5d_left", "Fig5d_right"
# sheet column -> internal name
LEFT_COLUMNS = {"sc_lr_score": "sc_raw", "visium_median_jaccard": "vis_raw",
                "x_sc_scaled": "sc_strength", "y_visium_scaled": "vis_strength"}
GROUP_PAIR = "Endothelial-Mural"   # scatter hue; the sheet has no group_pair column, all pairs are Endothelial-Mural

# "Outputs" ----
PLOTS_DIR, TABLES_DIR = cfg.output_dirs("Fig5")
OUT_FIG = PLOTS_DIR / "sc_vs_visium_bothsig_scatter_by_cancer_q95scaled_platformlists.pdf"

TICK_FONTSIZE = 15
LABEL_FONTSIZE = 17
AXIS_LABEL_FONTSIZE = 21   # shared x / y axis labels
TITLE_FONTSIZE = 16
N_ROWS, N_COLS = 4, 4
MAIN_COLS = 3          # cancer panels fill columns 1-3 top to bottom (alphabetical)
EXTRA_SLOT = (0, 3)    # 13th cancer type: row 1, column 4
FOREST_ROWS = slice(1, 4)   # rows 2-4, column 4 (below the colorbar)
WSPACE, HSPACE = 0.08, 0.26  # gaps between panels, as fractions of panel width / height
MARGINS = dict(left=0.075, right=0.99, bottom=0.06, top=0.935)
COL4_RATIO = 1.45   # wider right column (13th panel / forest plot / colorbar)
PANEL_W, PANEL_H = 3.0, 3.0   # square cancer panels

KDE_LEVELS = np.linspace(0.05, 1.0, 5)
BAND_CMAP = plt.get_cmap("viridis", len(KDE_LEVELS) - 1)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Fig. 5d: SC-vs-Visium scatter of LR pairs significant in both platforms "
                    "(input: Source_Data/Source_Data_Fig5.xlsx, sheets 'Fig.5d left' / 'Fig.5d right').",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    return p.parse_args()


def main() -> None:
    parse_args()
    # every plotted point and the per-cancer statistics of the forest plot
    points = cfg.read_source(LEFT_KEY)
    per_cancer = cfg.read_source(RIGHT_KEY)
    both = points.rename(columns=LEFT_COLUMNS).copy()
    both["group_pair"] = GROUP_PAIR

    # both sheets as read
    points.to_csv(TABLES_DIR / "source_data_scatter_points.csv", index=False)
    per_cancer.to_csv(TABLES_DIR / "source_data_per_cancer_spearman.csv", index=False)

    sns.set_style("white")
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    # mathtext (exponents, italic P) in Arial as well
    plt.rcParams["mathtext.fontset"] = "custom"
    plt.rcParams["mathtext.rm"] = "Arial"
    plt.rcParams["mathtext.it"] = "Arial:italic"
    plt.rcParams["mathtext.bf"] = "Arial:bold"

    pair_order = sorted(both["group_pair"].dropna().astype(str).unique().tolist())
    pair_palette = dict(zip(pair_order, sns.color_palette("husl", n_colors=len(pair_order))))
    cancers = sorted(both["cancer_type"].astype(str).unique().tolist())
    if len(cancers) != N_ROWS * MAIN_COLS + 1:
        raise ValueError(f"{len(cancers)} cancer types but layout has {N_ROWS * MAIN_COLS + 1} panels")


    def panel_scatter(data: pd.DataFrame, ax: plt.Axes) -> None:
        can_kde = len(data) >= 15 and data["sc_strength"].nunique() > 1 and data["vis_strength"].nunique() > 1
        kde_kws = dict(data=data, x="sc_strength", y="vis_strength", levels=KDE_LEVELS, bw_adjust=0.9,
                       clip=((0, 1), (0, 1)), warn_singular=False, ax=ax)
        if can_kde:
            sns.kdeplot(fill=True, cmap="viridis", **kde_kws)
            cs = ax.collections[-1]
            cs.set_norm(BoundaryNorm(cs.levels, BAND_CMAP.N))
            cs.set_cmap(BAND_CMAP)
        sns.scatterplot(data=data, x="sc_strength", y="vis_strength", hue="group_pair", hue_order=pair_order,
                        palette=pair_palette, s=12, alpha=0.45, linewidth=0, ax=ax, legend=False)
        if can_kde:
            sns.kdeplot(fill=False, color="#222222", linewidths=1.0, **kde_kws)
        ax.grid(False)


    n_rows, n_cols = N_ROWS, N_COLS
    # figure width chosen so that every grid cell is square (then set_box_aspect leaves no slack
    # and the horizontal gap between panels is exactly WSPACE x panel width)
    fig_h = PANEL_H * n_rows
    cell_h = fig_h * (MARGINS["top"] - MARGINS["bottom"]) / (n_rows + (n_rows - 1) * HSPACE)
    fig_w = cell_h * ((n_cols - 1) + COL4_RATIO + (n_cols - 1) * WSPACE) / (MARGINS["right"] - MARGINS["left"])
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(n_rows, n_cols, width_ratios=[1] * (n_cols - 1) + [COL4_RATIO])
    axes: dict[str, plt.Axes] = {}
    ref_ax = None
    check: dict[str, tuple[float, int, int]] = {}   # self-check only: rho recomputed from the raw values, n, n above y = x
    placements = [(ct, i // MAIN_COLS, i % MAIN_COLS) for i, ct in enumerate(cancers[: N_ROWS * MAIN_COLS])]
    placements.append((cancers[-1], *EXTRA_SLOT))
    for ct, r, c in placements:
        ax = fig.add_subplot(gs[r, c], sharex=ref_ax, sharey=ref_ax)
        ref_ax = ref_ax or ax
        d = both[both["cancer_type"] == ct]
        panel_scatter(d, ax)
        ax.plot([0, 1], [0, 1], color="white", linestyle="--", linewidth=1.2, zorder=5,   # y = x, white with a dark halo
                path_effects=[pe.Stroke(linewidth=2.4, foreground="black"), pe.Normal()])
        rho, _ = spearmanr(d["sc_raw"], d["vis_raw"])
        n_above = int((d["vis_strength"] > d["sc_strength"]).sum())   # points above y = x (scaled values)
        check[ct] = (rho, len(d), n_above)
        ax.set_title(f"{ct} (n = {len(d)})", size=TITLE_FONTSIZE)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="both", labelsize=TICK_FONTSIZE)
        ax.set_box_aspect(1)
        ax.label_outer()
        if (r, c) == EXTRA_SLOT:          # wider column: anchor to the left so the gap matches the other panels
            ax.set_anchor("W")
            ax.tick_params(axis="x", labelbottom=False)
        axes[ct] = ax
    ref_ax.set_xlim(-0.03, 1.03)
    ref_ax.set_ylim(-0.03, 1.03)
    ref_ax.set_xticks([0, 0.5, 1])
    ref_ax.set_yticks([0, 0.5, 1])
    ref_ax.set_xticklabels(["0", "0.5", "1"])   # short labels: no collision between neighbouring panels
    ref_ax.set_yticklabels(["0", "0.5", "1"])

    fig.supxlabel("scRNA-seq LR score (scaled)", fontsize=AXIS_LABEL_FONTSIZE, y=0.008)
    fig.supylabel("Visium LR co-localization (scaled)", fontsize=AXIS_LABEL_FONTSIZE, x=0.004)
    # ---- right column: forest plot (rows 2-4) and density colorbar (beside the row-1 panel) ----
    ax_f = fig.add_subplot(gs[FOREST_ROWS, n_cols - 1])   # rows 2-4 of the right column
    # rho, 95% CI (Fisher z) and n from the 'Fig.5d right' sheet; rows ordered by rho ascending (bottom to top)
    sd = per_cancer.sort_values("spearman_rho", kind="stable").reset_index(drop=True)
    order = sd["cancer_type"].tolist()
    rhos = sd["spearman_rho"].to_numpy(dtype=float)
    ns = sd["n_pairs"].to_numpy(dtype=int)
    lo, hi = sd["ci95_low"].to_numpy(dtype=float), sd["ci95_high"].to_numpy(dtype=float)
    n_above = sd["n_above_identity"].to_numpy(dtype=int)
    # self-check against the Source Data (not used for drawing)
    rho_diff = max(abs(check[c][0] - r) for c, r in zip(order, rhos))
    n_diff = max(max(abs(check[c][1] - n), abs(check[c][2] - k)) for c, n, k in zip(order, ns, n_above))
    print(f"self-check vs 'Fig.5d right': max |rho recomputed from sc_lr_score / visium_median_jaccard - spearman_rho| "
          f"= {rho_diff:.3e}; max |n_pairs / n_above_identity difference| = {n_diff}")
    print(f"pairs above y = x: pooled {n_above.sum() / ns.sum():.3f} ({n_above.sum()}/{ns.sum()}); "
          f"per cancer {sd.frac_above_identity.min():.2f}-{sd.frac_above_identity.max():.2f}")
    y = np.arange(len(order))
    for xv in (0, 0.25, 0.5):   # dotted guides
        ax_f.axvline(xv, color="#999999", linestyle=":", linewidth=0.8)
    ax_f.errorbar(rhos, y, xerr=[rhos - lo, hi - rhos], fmt="o", color="#d62728", ecolor="#444444",
                  elinewidth=0.9, capsize=2.5, markersize=6.5)
    ax_f.set_yticks(y)
    ax_f.set_yticklabels(order, fontsize=TICK_FONTSIZE)
    ax_f.set_xlabel("Spearman ρ (95% CI)", fontsize=LABEL_FONTSIZE)
    ax_f.set_title("Per-cancer correlation", size=TITLE_FONTSIZE)
    ax_f.tick_params(axis="x", labelsize=TICK_FONTSIZE)
    ax_f.set_xlim(min(-0.12, lo.min() - 0.06), max(0.7, hi.max() + 0.08))


    ax_f.set_ylim(-0.6, len(order) - 1 + 0.6)
    ax_f.set_xticks([0, 0.25, 0.5])
    for sp in ("top", "right"):
        ax_f.spines[sp].set_visible(False)
    ax_f.grid(False)


    fig.suptitle("Endothelial-Mural LR pairs significant in both platforms", fontsize=TITLE_FONTSIZE, y=0.985)
    # fixed, small gaps between panels (no tight_layout: inner panels carry no tick labels)
    fig.subplots_adjust(wspace=WSPACE, hspace=HSPACE, **MARGINS)
    # forest plot: shifted right within rows 2-4 to leave room for the cancer names
    pf = ax_f.get_position()
    ax_f.set_position([pf.x0 + 0.22 * pf.width, pf.y0, 0.78 * pf.width, pf.height])   # full rows 2-4; gap to the scatter column

    # vertical colorbar to the right of the 13th panel, matching its height
    pu = axes[cancers[-1]].get_position()
    fig_w_in = fig.get_size_inches()[0]
    cax = fig.add_axes([pu.x1 + 0.07 / fig_w_in, pu.y0, 0.14 / fig_w_in, pu.height])
    cbar = fig.colorbar(ScalarMappable(norm=BoundaryNorm(KDE_LEVELS, BAND_CMAP.N), cmap=BAND_CMAP), cax=cax, ticks=KDE_LEVELS)
    cbar.set_ticklabels([f"{(1 - p) * 100:.0f}%" for p in KDE_LEVELS])
    cbar.ax.tick_params(labelsize=TICK_FONTSIZE)
    cbar.set_label("Density", fontsize=LABEL_FONTSIZE)   # rotated label to the right of the colorbar

    fig.savefig(OUT_FIG, dpi=700, bbox_inches="tight")
    plt.close(fig)
    print("saved:", OUT_FIG)
    print("saved:", TABLES_DIR / "source_data_per_cancer_spearman.csv")
    print("saved:", TABLES_DIR / "source_data_scatter_points.csv")


if __name__ == "__main__":
    main()
