#!/usr/bin/env python3
"""Fig. 4a/b: per-slide median EMT signature scores (Epithelial, Mesenchymal, pEMT) by
compartment, 1 x 3 boxplots.

Input : Source_Data_Fig4.xlsx sheets 'Fig.4a' and 'Fig.4b' (cfg.read_source "Fig4a", "Fig4b"),
        268 slides x 3 compartments.
Output: Plots/fig4ab_signature_scores.pdf; Tables/fig4ab_statistics.csv, source_data_fig4a.csv,
        source_data_fig4b.csv under Output/Figure4/.
Compartments are subsets of the same section, so the on-plot p values are paired two-sided Wilcoxon
signed-rank tests on slide-matched medians, plus a Friedman test per signature; both are recomputed
from the sheet values (scipy).
Plotting code is that of Analysis/Run_Fig4_Statistics/01_Fig4ab_signature_panels.py, so the PDF matches.
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
from scipy.stats import friedmanchisquare, wilcoxon

plt.rcParams.update({
    "font.family": "Arial", "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})

COMP_COLORS = cfg.COMP_COLORS
ORDER = cfg.COMP_ORDER
LABELS = cfg.COMP_LABELS
CODES = {v: k for k, v in LABELS.items()}          # sheet column -> compartment code
SOURCE_SHEETS = {"Epithelial": "Fig4a", "Mesenchymal": "Fig4a", "pEMT": "Fig4b"}   # signature -> cfg.SHEETS key


def fmtp(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 1e-3 else f"{p:.3f}"


def facet(ax, m: pd.DataFrame, title: str, stats_rows: list, panel: str) -> None:
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
    sub = m[ORDER].dropna()
    try:
        p_fr = float(friedmanchisquare(*[sub[c] for c in ORDER]).pvalue)
    except ValueError:
        p_fr = 1.0
    stats_rows.append({"panel": panel, "signature": title, "comparison": "Mal-Bdy-Normal",
                       "test": "friedman", "n_slides": len(sub), "pvalue": p_fr})
    vals = m[ORDER].to_numpy()
    vmax = float(np.nanquantile(vals, 0.999))
    span = vmax - float(np.nanquantile(vals, 0.001)) or 1.0
    for (a, b), h in [((0, 1), 0.06), ((1, 2), 0.16), ((0, 2), 0.26)]:
        pair = m[[ORDER[a], ORDER[b]]].dropna()
        diff = pair.iloc[:, 0] - pair.iloc[:, 1]
        p = 1.0 if np.all(diff == 0) else float(wilcoxon(pair.iloc[:, 0], pair.iloc[:, 1]).pvalue)
        stats_rows.append({"panel": panel, "signature": title,
                           "comparison": f"{ORDER[a]}-{ORDER[b]}",
                           "test": "wilcoxon_signed_rank", "n_slides": len(pair), "pvalue": p})
        y = vmax + span * h
        ax.plot([a, a, b, b], [y, y + span * 0.012, y + span * 0.012, y], c="#333", lw=0.8)
        ax.text(1.0, y + span * 0.022, f"$p$={fmtp(p)}", ha="center", fontsize=10.5)
    ax.set_xticks(range(3))
    ax.set_xticklabels([LABELS[c] for c in ORDER], fontsize=15, rotation=45, ha="right")
    ax.tick_params(axis="y", labelsize=15)
    ax.set_title(title, fontsize=16)


def slide_table(sheets: dict[str, pd.DataFrame], name: str) -> pd.DataFrame:
    """Sheet -> slide x compartment table (index sample_id, columns ORDER codes), sheet row order kept."""
    df = sheets[SOURCE_SHEETS[name]]
    m = (df.loc[df["signature"] == name]
           .set_index("sample_id")[[LABELS[c] for c in ORDER]]
           .rename(columns=CODES))
    return m.dropna(subset=ORDER)


def main() -> None:
    # "Inputs" ----
    ap = argparse.ArgumentParser(description="Fig 4a/b paired signature-score panels from Source Data.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.parse_args()

    sheets = {key: cfg.read_source(key) for key in ("Fig4a", "Fig4b")}   # Source_Data_Fig4.xlsx 'Fig.4a' / 'Fig.4b'
    for key, df in sheets.items():
        assert df["sample_id"].nunique() == 268, (key, df["sample_id"].nunique())
    print(f"Source Data: Fig.4a {len(sheets['Fig4a'])} rows "
          f"({', '.join(sheets['Fig4a']['signature'].unique())}) | Fig.4b {len(sheets['Fig4b'])} rows "
          f"({', '.join(sheets['Fig4b']['signature'].unique())}); 268 slides x 3 compartments")

    # "Outputs" ----
    plots_dir, tables_dir = cfg.output_dirs("Fig4")   # Output/Figure4/{Plots,Tables}

    stats_rows = []
    # ---- Fig 4a/b combined: Epithelial + Mesenchymal + pEMT (equal facet boxes) ----
    fig, axes = plt.subplots(1, 3, figsize=(8.6, 4.4))
    src_a = []
    for ax, name in zip(axes, ("Epithelial", "Mesenchymal", "pEMT")):
        panel = "fig4b" if name == "pEMT" else "fig4a"
        m = slide_table(sheets, name)
        facet(ax, m, name, stats_rows, panel)
        s = m[ORDER].rename(columns=LABELS).reset_index()
        s.insert(0, "signature", name)
        if panel == "fig4a":
            src_a.append(s)
        else:
            s.to_csv(tables_dir / "source_data_fig4b.csv", index=False)
    axes[0].set_ylabel("Median score", fontsize=17)
    fig.tight_layout()
    fig.savefig(plots_dir / "fig4ab_signature_scores.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    pd.concat(src_a).to_csv(tables_dir / "source_data_fig4a.csv", index=False)

    pd.DataFrame(stats_rows).to_csv(tables_dir / "fig4ab_statistics.csv", index=False)
    print("saved:", plots_dir, "|", tables_dir)


if __name__ == "__main__":
    main()
