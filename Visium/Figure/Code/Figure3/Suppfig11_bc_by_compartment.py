#!/usr/bin/env python3
"""Supplementary Fig 11: subtype betweenness centrality in the Malignant, Boundary and Normal co-enrichment
networks, grouped horizontal bars per major lineage.

Input : Source_Data_Supp_Fig.xlsx sheet Sup fig.11 (cfg.read_source key SuppFig11); one row per subtype
        (subtype, lineage, Malignant, Boundary, Normal).
Output: Output/SupplementaryFigure11/Plots/suppfig11_betweenness_by_compartment.pdf,
        Tables/suppfig11_betweenness_by_compartment.csv, Tables/source_data_suppfig11.csv.
Centralities are read, not recomputed; subtypes with no edge in a network have 0 there (no bar). Plotting code
is that of Analysis/Run_Fig3_CoEnrichment/06_Suppfig11_bc_by_compartment.py, so the PDF matches.
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

cfg.setup_matplotlib()

# "Outputs" ----
PLOTS, TABLES = cfg.output_dirs("SuppFig11")

COMP_COLORS = cfg.COMP_COLORS
COMP_LABELS = cfg.COMP_LABELS
ORDER = cfg.COMP_ORDER
LINEAGE_ORDER = cfg.LINEAGE_ORDER


def main() -> None:
    ap = argparse.ArgumentParser(description="Supp Fig 11 (BC across compartments) from Source Data.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.parse_args()

    src = cfg.read_source("SuppFig11")                       # subtype, lineage, Malignant, Boundary, Normal
    out = src[["subtype", "lineage"] + [COMP_LABELS[c] for c in ORDER]]
    out.to_csv(TABLES / "suppfig11_betweenness_by_compartment.csv", index=False)
    out.to_csv(TABLES / "source_data_suppfig11.csv", index=False)

    # subtype x compartment table with code columns ORDER (Mal/Bdy/Normal) + lineage, as in the Analysis script
    bc = out.set_index("subtype")[[COMP_LABELS[c] for c in ORDER]].copy()
    bc.columns = ORDER
    bc["lineage"] = out["lineage"].to_numpy()

    fig, axes = plt.subplots(3, 3, figsize=(16.5, 17))
    for k, lin in enumerate(LINEAGE_ORDER):
        ax = axes[k // 3][k % 3]
        sub = bc[bc.lineage == lin].copy()
        sub["mx"] = sub[ORDER].max(axis=1)
        sub = sub.sort_values("mx", ascending=False)     # largest at the bottom row
        y = range(len(sub))
        for off, comp in zip((0.27, 0.0, -0.27), ORDER):  # top-to-bottom: Mal, Bdy, Normal
            ax.barh([yy + off for yy in y], sub[comp], height=0.25,
                    color=COMP_COLORS[comp], edgecolor="#222", linewidth=0.4)
        ax.set_yticks(list(y))
        ax.set_yticklabels(sub.index, fontsize=10)
        ax.set_title(lin, fontsize=15)
        ax.set_xlabel("Betweenness Centrality", fontsize=12)
        ax.tick_params(axis="x", labelsize=10)
        ax.set_ylim(-0.7, len(sub) - 0.3)
        ax.set_axisbelow(True)
        ax.grid(axis="x", color="#E3E3E3", linewidth=0.7)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=COMP_COLORS[c], edgecolor="#222",
                             label=COMP_LABELS[c]) for c in ORDER]
    axes[2][1].axis("off")
    axes[2][2].axis("off")
    axes[2][1].legend(handles=handles, loc="center", frameon=False, fontsize=14)
    fig.tight_layout()
    fig.savefig(PLOTS / "suppfig11_betweenness_by_compartment.pdf",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("saved:", PLOTS / "suppfig11_betweenness_by_compartment.pdf")


if __name__ == "__main__":
    main()
