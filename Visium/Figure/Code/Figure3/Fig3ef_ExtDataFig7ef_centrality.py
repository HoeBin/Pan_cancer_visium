#!/usr/bin/env python3
"""Fig 3e/3f and Extended Data Fig 7e/7f: betweenness centrality, lineage median bars
(e), ranked subtype scatter (f).

Input : Source_Data_Fig3.xlsx sheets Fig.3e, Fig.3f (cfg.read_source keys Fig3e,
        Fig3f); Source_Data_Extended_Fig7.xlsx sheets ED Fig.7e, ED Fig.7f (ExtFig7e,
        ExtFig7f; BRCA-downsampled cohort).
Output: Plots/<fig3ef|extfig7ef>_centrality.pdf, Tables/source_data_<fig3e|fig3f|extfig7e|extfig7f>.csv under
        Output/Figure3 or Output/ExtendedFigure7.
Betweenness is not recomputed: panel e reindexes the sheet to cfg.LINEAGE_ORDER, panel f uses the rows in
sheet order (= rank) and labels the top --top-label nodes. Plotting code is that of
Analysis/Run_Fig3_CoEnrichment/02_Fig3ef_centrality.py, so the PDFs match. --panel fig3ef|extfig7ef|all.
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
import pandas as pd

plt.rcParams.update({
    "font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})

LINEAGE_COLORS = cfg.LINEAGE_COLORS
LINEAGE_ORDER = cfg.LINEAGE_ORDER

# --panel -> (panel-e sheet key, panel-f sheet key, output figure, e/f table stems, PDF name)
PANELS = {
    "fig3ef": ("Fig3e", "Fig3f", "Fig3", "fig3e", "fig3f", "fig3ef_centrality.pdf"),
    "extfig7ef": ("ExtFig7e", "ExtFig7f", "ExtFig7", "extfig7e", "extfig7f", "extfig7ef_centrality.pdf"),
}


def draw_centrality(e_key: str, f_key: str, figure: str, stem_e: str, stem_f: str, pdf_name: str,
                    top_label: int) -> None:
    """One 1x2 panel pair: panel e from the median sheet, panel f from the ranked sheet."""
    figdir, tabdir = cfg.output_dirs(figure)
    med_sheet = cfg.read_source(e_key)        # lineage, median_betweenness
    tab = cfg.read_source(f_key)              # rank, node, betweenness, lineage (sheet order = rank order)
    print(f"[{e_key}/{f_key}] lineages={len(med_sheet)}  nodes={len(tab)}")

    # Series indexed by lineage in LINEAGE_ORDER, as in the Analysis script
    med = med_sheet.set_index("lineage")["median_betweenness"].reindex(LINEAGE_ORDER)
    pd.DataFrame({"lineage": med.index, "median_betweenness": med.values}
                 ).to_csv(tabdir / f"source_data_{stem_e}.csv", index=False)
    tab.to_csv(tabdir / f"source_data_{stem_f}.csv", index=False)

    # ---- Fig 3e + 3f combined (1x2, identical axes geometry) ----
    fig, (axe, axf) = plt.subplots(1, 2, figsize=(11.8, 4.8))

    axe.bar(range(len(med)), med.values, color=[LINEAGE_COLORS[l] for l in med.index],
            edgecolor="#333", linewidth=0.6)
    axe.set_xticks(range(len(med)))
    axe.set_xticklabels(med.index, rotation=40, ha="right", fontsize=17)
    axe.set_ylabel("Median Betweenness Centrality", fontsize=16)
    axe.tick_params(axis="y", labelsize=14)
    for s_ in ("top", "right"):
        axe.spines[s_].set_visible(False)

    axf.scatter(range(1, len(tab) + 1), tab.betweenness, s=42,
                c=[LINEAGE_COLORS[l] for l in tab.lineage], edgecolors="#333", linewidths=0.5)
    n_lab = min(top_label, len(tab))
    if n_lab:
        ymax = float(tab.betweenness.iloc[0])
        lx = 0.32 * len(tab)
        for i in range(n_lab):
            axf.annotate(tab.node.iloc[i], xy=(i + 1, float(tab.betweenness.iloc[i])),
                         xytext=(lx, ymax * (1.0 - 0.075 * i)), va="center", fontsize=13.5,
                         arrowprops={"arrowstyle": "-", "lw": 0.6, "color": "#999",
                                     "shrinkA": 2, "shrinkB": 3})
    axf.set_xlabel("Rank", fontsize=19)
    axf.set_ylabel("Betweenness Centrality", fontsize=16)
    axf.set_xticks([])
    axf.tick_params(axis="y", labelsize=14)
    handles = [plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=c,
                          markeredgecolor="#333", markersize=8, label=l)
               for l, c in LINEAGE_COLORS.items()]
    axf.legend(handles=handles, frameon=False, fontsize=12.5, loc="upper right")
    for s_ in ("top", "right"):
        axf.spines[s_].set_visible(False)
    fig.tight_layout()
    fig.savefig(figdir / pdf_name, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("top 10 betweenness:")
    print(tab.head(10).round(4).to_string(index=False))
    print("saved:", figdir / pdf_name)


def main() -> None:
    ap = argparse.ArgumentParser(description="Fig 3e/3f and Extended Data Fig 7e/7f betweenness plots "
                                             "from the Source Data sheets.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--panel", choices=("fig3ef", "extfig7ef", "all"), default="all",
                    help="which panel pair to draw (fig3ef: full cohort; extfig7ef: BRCA-downsampled)")
    ap.add_argument("--top-label", type=int, default=10, help="number of top-ranked nodes to label in panel f")
    args = ap.parse_args()

    panels = list(PANELS) if args.panel == "all" else [args.panel]
    for p in panels:
        e_key, f_key, figure, stem_e, stem_f, pdf_name = PANELS[p]
        draw_centrality(e_key, f_key, figure, stem_e, stem_f, pdf_name, args.top_label)


if __name__ == "__main__":
    main()
