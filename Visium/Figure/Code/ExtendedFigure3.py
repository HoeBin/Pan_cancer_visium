#!/usr/bin/env python3
"""Extended Data Fig 3: row z-score heatmaps of the median per-slide Pearson correlation between marker gene
expression and deconvolution abundance, one panel per lineage (a: 7 major lineages, b-h: subtypes).

Input : Source_Data_Extended_Fig3.xlsx sheet ED Fig.3 (cfg.read_source key ExtFig3);
        ~354k rows, about 1 min to read.
Output: Output/ExtendedFigure3/Plots/heatmap_corr_<Major|lineage>_median_pearson_zscore.pdf; Tables/
        source_data_ExtDataFig3.csv, extfig3_median_pearson.csv, extfig3_major_median_pearson.csv,
        extfig3_plotted_zscores.csv.
Recomputed from the sheet: median over slides per (gene, celltype) and the row z-score per gene across
columns (fixed +/-3 RdBu_r scale); cells where the row gene is the designated marker of the column are
outlined (cfg.EXTFIG3_PANELS).
"""
from __future__ import annotations

# "Libraries and paths" ----
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))   # Code/ for `common`
from common import config as cfg

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

plt.rcParams["font.family"] = "Arial"

FS_TITLE = 18       # figure title (global_anno)
FS_LABEL = 15       # gene / sub_anno tick labels
FS_TICK = 11        # colorbar tick numbers
FS_CBLABEL = 14     # colorbar label
CELL = 0.45         # cell size (inches)
VMAX_Z = 3.0        # fixed z-score scale; |z| > 3 saturates

PANELS = cfg.EXTFIG3_PANELS
MAJORS = cfg.EXTFIG3_MAJORS
PANEL_LETTER = cfg.EXTFIG3_PANEL_LETTER

# "Inputs" ----
SOURCE_KEY = "ExtFig3"          # Source_Data_Extended_Fig3.xlsx / 'ED Fig.3'

# "Outputs" ----
PLOTS, TABLES = cfg.output_dirs("ExtFig3")


def load_persample(args) -> pd.DataFrame:
    """slide-level Pearson r table (panel, global_anno, sample_id, gene, celltype, pearson_r) from the sheet."""
    ps = cfg.read_source(SOURCE_KEY)
    ps["pearson_r"] = pd.to_numeric(ps["pearson_r"])
    fixed = ps["celltype"].map(cfg.canonical_subtype)          # map sheet spellings to the config names
    changed = sorted(set(ps.loc[fixed != ps["celltype"], "celltype"]))
    if changed:
        print(f"[note] celltype spelling normalised via cfg.SUBTYPE_ALIASES: {changed}")
    ps["celltype"] = fixed
    return ps


def panel_specs():
    """(letter, title, genes, cols, marks, designated) per panel, in figure order a-h."""
    specs = []
    majors = [m for m, _ in MAJORS]
    genes = [g for _, g in MAJORS]
    specs.append(("a", "Major celltypes", genes, majors, [(i, i) for i in range(len(majors))], None))
    for gname, rows in PANELS:
        subs = [s for s, _, _ in rows]
        designated = {s: g for s, _, g in rows}
        genes = list(dict.fromkeys(g for _, _, g in rows))   # unique, in order
        marks = [(i, j) for i, g in enumerate(genes) for j, s in enumerate(subs) if designated[s] == g]
        specs.append((PANEL_LETTER[gname], gname, genes, subs, marks, designated))
    return sorted(specs, key=lambda t: t[0])


def plot_heatmap(gname, genes, subs, M, marks, vmax, out_pdf, png_copy=None,
                 cbar_label="Row z-score of median Pearson r"):
    """M: len(genes) x len(subs); marks: (row, col) cells to outline."""
    n_r, n_c = len(genes), len(subs)
    left_in = 0.075 * (FS_LABEL / 10) * max(len(g) for g in genes) + 0.35
    bottom_in = 0.058 * (FS_LABEL / 10) * max(len(s) for s in subs) + 0.55
    cbar_in = 1.05
    fig_w = CELL * n_c + left_in + cbar_in + 0.35
    fig_h = CELL * n_r + bottom_in + 0.55
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([left_in / fig_w, bottom_in / fig_h,
                       CELL * n_c / fig_w, CELL * n_r / fig_h])
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="equal")
    for i, j in marks:          # outline the designated marker cells
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                               edgecolor="black", lw=1.4))
    ax.set_xticks(range(n_c))
    ax.set_xticklabels(subs, rotation=45, ha="right", rotation_mode="anchor",
                       fontsize=FS_LABEL)
    ax.set_yticks(range(n_r))
    ax.set_yticklabels(genes, fontsize=FS_LABEL)
    ax.tick_params(length=2)
    ax.set_xticks(np.arange(-0.5, n_c), minor=True)
    ax.set_yticks(np.arange(-0.5, n_r), minor=True)
    ax.grid(which="minor", color="white", lw=0.6)
    ax.tick_params(which="minor", length=0)

    cax = fig.add_axes([(left_in + CELL * n_c + 0.25) / fig_w,
                        bottom_in / fig_h, 0.18 / fig_w, CELL * n_r / fig_h])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=FS_TICK, length=2)
    cb.set_label(cbar_label, fontsize=FS_CBLABEL)
    ax.set_title(gname, fontsize=FS_TITLE, pad=12)
    fig.savefig(out_pdf)
    if png_copy:
        fig.savefig(png_copy, dpi=110)
    plt.close(fig)
    print(f"saved {out_pdf}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Extended Data Fig. 3 heatmaps (median per-sample Pearson r, row z-score).",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    args = ap.parse_args()

    ps = load_persample(args)
    ps.to_csv(TABLES / "source_data_ExtDataFig3.csv", index=False)
    med_all = ps.groupby(["panel", "gene", "celltype"], sort=False)["pearson_r"].median()

    med_rows, major_rows, src_rows = [], [], []
    for letter, gname, genes, cols, marks, designated in panel_specs():
        M = np.full((len(genes), len(cols)), np.nan)
        for i, g in enumerate(genes):
            for j, s in enumerate(cols):
                M[i, j] = med_all.get((letter, g, s), np.nan)
                if letter == "a":
                    major_rows.append({"gene": g, "marker_of": MAJORS[i][0], "major": s, "median_pearson_r": M[i, j]})
                else:
                    med_rows.append({"global_anno": gname, "gene": g, "sub_anno": s,
                                     "median_pearson_r": M[i, j], "is_designated": designated[s] == g})
        # row-wise z-score (per gene across columns)
        mu = np.nanmean(M, axis=1, keepdims=True)
        sd = np.nanstd(M, axis=1, keepdims=True)
        sd[sd == 0] = np.nan
        Z = (M - mu) / sd
        for i, g in enumerate(genes):
            for j, s in enumerate(cols):
                src_rows.append({"panel": letter, "global_anno": gname, "gene": g, "celltype": s,
                                 "median_pearson_r": M[i, j], "row_zscore": Z[i, j],
                                 "is_designated": (i, j) in marks})
        safe = "Major" if letter == "a" else gname.replace(" ", "_").replace("/", "_")
        plot_heatmap(gname, genes, cols, Z, marks, VMAX_Z,
                     PLOTS / f"heatmap_corr_{safe}_median_pearson_zscore.pdf")

    # median table in EXTFIG3_PANELS order
    med = pd.DataFrame(med_rows)
    order = {g: i for i, (g, _) in enumerate(PANELS)}
    med = med.iloc[np.argsort(med["global_anno"].map(order).values, kind="stable")].reset_index(drop=True)
    med.to_csv(TABLES / "extfig3_median_pearson.csv", index=False)
    pd.DataFrame(major_rows).to_csv(TABLES / "extfig3_major_median_pearson.csv", index=False)
    pd.DataFrame(src_rows).to_csv(TABLES / "extfig3_plotted_zscores.csv", index=False)
    print("saved tables:", TABLES)


if __name__ == "__main__":
    main()
