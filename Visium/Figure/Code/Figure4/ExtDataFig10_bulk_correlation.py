#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Ext. Data Fig. 10: rank dotplots of the bulk LR product vs HSPA6+ Macro x tCAF product Spearman
correlation per cancer type.

Input : Source_Data_Extended_Fig10.xlsx sheet 'ED Fig.10' (cfg.read_source
        "ExtFig10"), one row per plotted point.
Output: Plots/<cancer>.pdf (HNCA, KICA, LICA, LUCA, PACA); Tables/source_data_ExtDataFig10.csv under
        Output/ExtendedFigure10/ (--outdir changes the root; old PDFs in Plots/ are removed first).
Correlations, q values and ranks are the sheet values and highlighted pairs come from the Highlighted
column; nothing is recomputed, and correlation_summary.csv of the Analysis script is not written (the
sheet has no Spearman P column).
Needs adjustText. Plotting code is that of
Analysis/Run_Fig4_Statistics/04_ExtDataFig10_bulk_correlation.py, so the PDFs match.
"""
from __future__ import annotations

# "Libraries and paths" ----
import argparse
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # Code/ for `common`
from common import config as cfg

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("KMP_AFFINITY", "disabled")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib import font_manager
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

# "Inputs" ----
CELLTYPE_PAIR = "HSPA6+ Macro - tCAF"
PRIORITIZED_PAIRS = ["PLAU-PLAUR", "FN1-PLAUR", "SPP1-ITGA5", "SERPINE1-PLAUR"]
SIGNIFICANCE_COL = "q_value_bh_per_cancer"
# sheet column -> internal name
SOURCE_COLUMNS = {
    "Cancer type": "cancer_group", "LR pair": "lr_pair", "Spearman correlation": "correlation",
    "No. of tumors": "n", "BH-adjusted q value (FDR)": SIGNIFICANCE_COL,
}

# "Outputs" ----
# Default root is Output/ExtendedFigure10/ (override with --outdir):
#   <outdir>/Plots/<cancer>.pdf ; <outdir>/Tables/source_data_ExtDataFig10.csv
PLOTS_DIR, TABLES_DIR = cfg.output_dirs("ExtFig10")
DEFAULT_OUTDIR = PLOTS_DIR.parent

# Panel typography
FS_PAIR_LABEL = 10
FS_AXIS_LABEL = 16
FS_TICK = 13
FS_TITLE = 14
FS_CBAR_LABEL = 13
FS_CBAR_TICK = 12


def configure_plot_style() -> None:
    try:
        font_manager.findfont("Arial", fallback_to_default=False)
        family = "Arial"
    except ValueError:
        family = plt.rcParams["font.sans-serif"][0]
        print(f"[configure_plot_style] Arial not found; using '{family}' instead.", flush=True)
    plt.rcParams.update({
        "font.family": family, "font.sans-serif": [family],
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    if family == "Arial":
        plt.rcParams.update({
            "mathtext.fontset": "custom", "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
        })


def sanitize_name(name: str) -> str:
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(name))


def highlight_map_from_source(src: pd.DataFrame) -> dict[str, list[str]]:
    """{cancer type: [prioritized pairs highlighted in that cancer type]} from the sheet's Highlighted column."""
    highlighted = src["Highlighted"].map(lambda v: str(v).strip().lower() == "true")
    hl = src[highlighted]
    cancer_to_pairs: dict[str, list[str]] = {}
    for pair in PRIORITIZED_PAIRS:
        for cancer in hl.loc[hl["LR pair"] == pair, "Cancer type"]:
            cancer_to_pairs.setdefault(cancer, [])
            if pair not in cancer_to_pairs[cancer]:
                cancer_to_pairs[cancer].append(pair)
    other = sorted(set(hl["LR pair"]) - set(PRIORITIZED_PAIRS))
    if other:
        print(f"  ! highlighted pairs outside PRIORITIZED_PAIRS are ignored: {other}", flush=True)
    return cancer_to_pairs


def plot_cancer_pairs(res_df: pd.DataFrame, cancer: str, pairs: list[str], out_path: Path) -> None:
    """Rank dotplot of one cancer type, highlighting and labelling ``pairs``."""
    df = res_df.sort_values("correlation", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1

    fig, ax = plt.subplots(figsize=(6, 5))
    nl_p = -np.log10(df[SIGNIFICANCE_COL].clip(lower=1e-300).values)
    norm = Normalize(vmin=0.0, vmax=max(float(np.nanmax(nl_p)), 3.0))
    bg_scatter = ax.scatter(df["rank"], df["correlation"], c=nl_p, cmap="viridis", norm=norm,
                            s=10, edgecolors="none", alpha=0.85)
    ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.6, linestyle="--")

    selected = df[df["lr_pair"].isin(pairs)].copy()
    missing = [p for p in pairs if p not in set(selected["lr_pair"])]
    if missing:
        print(f"  ! {cancer}: pairs not present, skipped -> {missing}", flush=True)
    if not selected.empty:
        selected = selected.sort_values("correlation", ascending=False).reset_index(drop=True)
        selected_nl_p = -np.log10(selected[SIGNIFICANCE_COL].clip(lower=1e-300).values)
        sel_scatter = ax.scatter(selected["rank"], selected["correlation"],
                                 c=plt.cm.viridis(norm(selected_nl_p)),
                                 s=60, edgecolors="red", linewidths=1.4, zorder=5)
        ax.margins(x=0.12, y=0.26)
        ymax = float(df["correlation"].max())
        ymin = float(df["correlation"].min())
        yr = (ymax - ymin) or 1.0
        rank_max = float(df["rank"].max())
        n_sel = len(selected)
        label_x = rank_max * 0.48
        y_hi = ymax + 0.28 * yr
        y_lo = ymax - 0.18 * yr
        texts = []
        for i, row in selected.iterrows():
            ty = y_hi if n_sel == 1 else y_hi - (y_hi - y_lo) * i / (n_sel - 1)
            texts.append(ax.text(label_x, ty, f"{row['lr_pair']} (#{int(row['rank'])})",
                                 fontsize=FS_PAIR_LABEL, color="black", fontweight="bold",
                                 zorder=6, ha="left", va="center"))
        fig.canvas.draw()
        adjust_text(
            texts, ax=ax,
            target_x=selected["rank"].values, target_y=selected["correlation"].values,
            objects=[bg_scatter, sel_scatter],
            arrowprops=dict(arrowstyle="-", color="red", lw=0.7),
            expand=(1.6, 1.9), force_text=(0.6, 1.1), force_static=(0.5, 0.9),
            force_pull=(0.002, 0.002), ensure_inside_axes=True, explode_radius=0, time_lim=12,
        )

    ax.set_xlabel("Rank (correlation, descending)", fontsize=FS_AXIS_LABEL)
    ax.set_ylabel("Spearman correlation\n(LR product x CT product)", fontsize=FS_AXIS_LABEL)
    ax.tick_params(axis="both", labelsize=FS_TICK)
    ax.grid(axis="y", alpha=0.25)
    ax.set_title(f"{cancer} (n={int(df['n'].iloc[0])})  |  {CELLTYPE_PAIR}", fontsize=FS_TITLE)
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap="viridis"), ax=ax, pad=0.01, fraction=0.04, shrink=0.85)
    cbar.set_label(r"-log10(adjusted $P$-value)", fontsize=FS_CBAR_LABEL)
    cbar.ax.tick_params(labelsize=FS_CBAR_TICK)
    fig.tight_layout()
    fig.savefig(out_path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Bulk LR-product x HSPA6+ Macro-tCAF correlation rank dotplots (Extended Data Fig. 10) "
                    "from the paper's Source Data.")
    ap.add_argument("--outdir", default=str(DEFAULT_OUTDIR),
                    help="output root: <cancer>.pdf under <outdir>/Plots/, source_data_ExtDataFig10.csv under "
                         "<outdir>/Tables/ (default: %(default)s)")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir).resolve()
    plots_dir, tables_dir = outdir / "Plots", outdir / "Tables"
    configure_plot_style()
    plots_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    for stale_pdf in plots_dir.glob("*.pdf"):
        stale_pdf.unlink()

    # ---- Source Data (one row per plotted point) ----
    src = cfg.read_source("ExtFig10")          # Source_Data_Extended_Fig10.xlsx 'ED Fig.10'
    missing = [c for c in [*SOURCE_COLUMNS, "Rank", "Highlighted"] if c not in src.columns]
    if missing:
        raise SystemExit(f"[ExtFig10] Source Data sheet lacks columns: {missing}")
    src.to_csv(tables_dir / "source_data_ExtDataFig10.csv", index=False)
    print(f"Source Data rows: {len(src)} -> {tables_dir / 'source_data_ExtDataFig10.csv'}")

    cancer_to_pairs = highlight_map_from_source(src)
    cancers = sorted(src["Cancer type"].unique())
    print("Cancer types with a prioritized pair associated with shorter survival:")
    for cancer in cancers:
        print(f"  {cancer}: {', '.join(cancer_to_pairs.get(cancer, [])) or '(none highlighted)'}")

    results = src.rename(columns=SOURCE_COLUMNS)
    for cancer in cancers:
        res = results[results["cancer_group"] == cancer].sort_values("Rank").reset_index(drop=True)
        plot_cancer_pairs(res, cancer, cancer_to_pairs.get(cancer, []), plots_dir / f"{sanitize_name(cancer)}.pdf")
        print(f"{cancer}: n={int(res['n'].iloc[0])} | LR pairs plotted={len(res)} | saved {cancer}.pdf")


if __name__ == "__main__":
    main()
