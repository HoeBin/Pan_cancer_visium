#!/usr/bin/env python3
"""Fig. 4d: LR pairs ranked by pan-cancer median Jaccard index with HSPA6+ Macro-tCAF in
the Boundary compartment.

Input : Source_Data_Fig4.xlsx sheet 'Fig.4d' (cfg.read_source "Fig4d"), one row per plotted point.
Output: Plots/Fig4d_cancertype.pdf (+ .png); Tables/source_data_Fig4d.csv under Output/Figure4/.
Median Jaccard, replication P, BH q and Rank are the sheet values; nothing is recomputed.
Of the top --n-top ranked pairs only those with replication q < --alpha are labelled.
Plotting code is that of Analysis/Run_Fig4_Statistics/02_Fig4d_cancertype_perm_test.py, so the PDF matches.
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

# "Inputs" ----
# sheet column -> internal name
SOURCE_COLUMNS = {
    "LR pair": "lr_pair",
    "Median Jaccard index": "median_jaccard",
    "No. of cancer types tested": "n_cancers_tested",
    "No. of cancer types with P < 0.01": "n_cancers_significant",
    "Replication P value": "replication_p",
    "BH-adjusted q value (FDR)": "replication_q_bh",
    "Rank": "rank",
}


def main() -> None:
    ap = argparse.ArgumentParser(description="Fig 4d ranked co-localization scatter from Source Data.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--alpha", type=float, default=0.05, help="replication q threshold for labeling")
    ap.add_argument("--n-top", type=int, default=10)
    args = ap.parse_args()

    # "Outputs" ----
    plots, tables = cfg.output_dirs("Fig4")   # Output/Figure4/{Plots,Tables}

    # ---- Source Data for Fig. 4d (one row per plotted point) ----
    src = cfg.read_source("Fig4d")            # Source_Data_Fig4.xlsx 'Fig.4d'
    missing = [c for c in SOURCE_COLUMNS if c not in src.columns]
    if missing:
        raise SystemExit(f"[Fig4d] Source Data sheet lacks columns: {missing}")
    src.to_csv(tables / "source_data_Fig4d.csv", index=False)
    full = src.rename(columns=SOURCE_COLUMNS)
    print(f"Source Data rows: {len(full)} | ranked pairs: {int(full['rank'].max())} "
          f"| q<{args.alpha}: {int((full['replication_q_bh'] < args.alpha).sum())}")

    top = full[full["rank"] <= args.n_top].copy()
    top["labeled"] = top["replication_q_bh"] < args.alpha
    top["vote"] = top["n_cancers_significant"].astype(int).astype(str) + "/" + top["n_cancers_tested"].astype(int).astype(str)
    print(top[["rank", "lr_pair", "median_jaccard", "vote", "replication_p",
               "replication_q_bh", "labeled"]].to_string(index=False))
    dropped = top.loc[~top["labeled"], "lr_pair"].tolist()
    if dropped:
        print(f"top-{args.n_top} pairs NOT labeled (q >= {args.alpha}): {dropped}")

    # ---- plot ----
    lab = top[top["labeled"]].sort_values("rank")
    fig, ax = plt.subplots(figsize=(6.0, 3.6), dpi=200)
    ax.scatter(full["rank"], full["median_jaccard"], s=10, color="black", zorder=2)
    ax.scatter(lab["rank"], lab["median_jaccard"], s=34, facecolor="black",
               edgecolor="red", linewidth=1.6, zorder=3)
    n_rank = int(full["rank"].max())
    x_text = n_rank * 0.52
    ymax = float(full["median_jaccard"].max())
    if len(lab):
        y_slots = np.linspace(ymax * 0.98, ymax * 0.30, num=len(lab))
        for (_, r), y_lab in zip(lab.iterrows(), y_slots):
            ax.plot([r["rank"], x_text - n_rank * 0.01], [r["median_jaccard"], y_lab],
                    color="0.35", lw=0.6, zorder=1)
            ax.text(x_text, y_lab, r["lr_pair"], fontsize=8.5, va="center", ha="left")
    pad = n_rank * 0.005                                   # axis ends at the last-ranked dot
    ax.set_xlim(1 - pad, n_rank + pad)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Rank", fontsize=12)
    ax.set_ylabel("Median Co-localization score", fontsize=12)
    ax.set_title("Boundary co-localization between\nHSPA6+ Macro–tCAF and LR pairs",
                 fontsize=12)
    ax.set_xticks([])
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(plots / "Fig4d_cancertype.pdf")
    fig.savefig(plots / "Fig4d_cancertype.png")
    print(f"saved {plots / 'Fig4d_cancertype.pdf'} (+.png), labeled={len(lab)}/{args.n_top}")


if __name__ == "__main__":
    main()
