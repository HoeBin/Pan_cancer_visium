#!/usr/bin/env python3
"""Fig 3c/3g/3h and Extended Data Fig 7c/7g/7h: per-slide median Jaccard by compartment (c) and gateway-anchor
partner Jaccard by compartment (g: Cycling_Bgc, h: HSPA6+_Macro).

Input : Source_Data_Fig3.xlsx sheets Fig.3c/3g/3h (cfg.read_source keys Fig3c, Fig3g,
        Fig3h); Source_Data_Extended_Fig7.xlsx sheets ED Fig.7c/7g/7h (ExtFig7c, ExtFig7g,
        ExtFig7h; BRCA-downsampled cohort).
Output: Plots/<prefix>c_per_sample_J.pdf, <prefix>g_Cycling_Bgc_partners.pdf,
        <prefix>h_HSPA6+_Macro_partners.pdf, <prefix>gh_legend.pdf; Tables/source_data_<prefix>{c,g,h}.csv,
        <prefix>cgh_statistics.csv (prefix fig3 | extfig7).
Only the paired two-sided Wilcoxon signed-rank tests on slide-matched values are recomputed; rows are used in
sheet order (partner order in g/h). Plotting and test code is that of
Analysis/Run_Fig3_CoEnrichment/03_Fig3cgh.py. --panel, --panels.
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
    "font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})

COMP_COLORS = cfg.COMP_COLORS
ORDER = cfg.COMP_ORDER
LABELS = cfg.COMP_LABELS
ANCHORS = {"g": "Cycling_Bgc", "h": "HSPA6+_Macro"}

# "Inputs / outputs per panel set" ----
# panel set -> output figure key, file-name prefix, Source Data keys per sub-panel
PANEL_SETS = {
    "fig3":    {"figure": "Fig3",    "prefix": "fig3",
                "sheets": {"c": "Fig3c", "g": "Fig3g", "h": "Fig3h"}},
    "extfig7": {"figure": "ExtFig7", "prefix": "extfig7",
                "sheets": {"c": "ExtFig7c", "g": "ExtFig7g", "h": "ExtFig7h"}},
}


STATS: list = []


def record(panel: str, group: str, comparison: str, test: str, n: int, p: float) -> None:
    STATS.append({"panel": panel, "group": group, "comparison": comparison,
                  "test": test, "n_slides": n, "pvalue": p})


def fmtp(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 1e-3 else f"{p:.3f}"


def paired_wilcoxon(piv: pd.DataFrame, a: str, b: str) -> tuple[float, int]:
    """Two-sided Wilcoxon signed-rank on slides having both compartments."""
    sub = piv[[a, b]].dropna()
    if len(sub) < 3:
        return np.nan, len(sub)
    diff = sub[a].to_numpy() - sub[b].to_numpy()
    if np.all(diff == 0):
        return 1.0, len(sub)
    return float(wilcoxon(sub[a], sub[b]).pvalue), len(sub)


def savefig(fig, outdir: Path, name: str, tight: bool = True) -> None:
    for ext in ("pdf",):
        fig.savefig(outdir / f"{name}.{ext}", dpi=300,
                    bbox_inches="tight" if tight else None)
    plt.close(fig)
    print("saved:", outdir / f"{name}.pdf")


def out_name(name: str, prefix: str) -> str:
    """File name for a panel: 'fig3c_per_sample_J' -> 'extfig7c_per_sample_J' when prefix == 'extfig7'."""
    return prefix + name[len("fig3"):]


def compartment_matrix(df: pd.DataFrame, index) -> pd.DataFrame:
    """Sheet columns Malignant/Boundary/Normal -> slide x compartment table with code columns ORDER (Mal/Bdy/Normal),
    as the Analysis script builds it from the merged table."""
    m = df.set_index(index)[[LABELS[c] for c in ORDER]].copy()
    m.columns = ORDER
    return m


def compartment_legend(outdir: Path, prefix: str = "fig3") -> None:
    """Standalone horizontal compartment legend (shared by Fig 3g/h)."""
    fig = plt.figure(figsize=(6.0, 0.6))
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=COMP_COLORS[c], edgecolor="#111",
                             label=LABELS[c]) for c in ORDER]
    fig.legend(handles=handles, frameon=False, fontsize=14, ncol=3, loc="center")
    name = out_name("fig3gh_legend", prefix)
    fig.savefig(outdir / f"{name}.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("saved:", outdir / f"{name}.pdf")


def fig3c(sheet: pd.DataFrame, outdir: Path, source_dir: Path, prefix: str = "fig3") -> None:
    name = "fig3c_per_sample_J"
    if "sample" in sheet.columns and "sample_id" not in sheet.columns:   # Fig.3c sheet header is `sample`
        sheet = sheet.rename(columns={"sample": "sample_id"})
    m = compartment_matrix(sheet, "sample_id").dropna(subset=ORDER)      # slides with all 3 compartments
    src = m[ORDER].rename(columns=LABELS).reset_index()
    src.to_csv(source_dir / f"source_data_{prefix}c.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    rng = np.random.default_rng(0)
    for i, c in enumerate(ORDER):
        v = m[c].to_numpy()
        ax.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v, s=7, c="#999", alpha=0.45, linewidths=0, zorder=1)
        ax.boxplot([v], positions=[i], widths=0.5, showfliers=False, patch_artist=True, zorder=2,
                   boxprops={"facecolor": "none", "edgecolor": COMP_COLORS[c], "linewidth": 1.7},
                   medianprops={"color": COMP_COLORS[c], "linewidth": 1.9},
                   whiskerprops={"color": COMP_COLORS[c]}, capprops={"color": COMP_COLORS[c]})
    ymax = float(np.nanquantile(m.to_numpy(), 0.999))
    for (a, b), h in [((0, 1), 1.03), ((1, 2), 1.13), ((0, 2), 1.23)]:
        p, n = paired_wilcoxon(m, ORDER[a], ORDER[b])
        record(name, "all_slides", f"{ORDER[a]}-{ORDER[b]}", "wilcoxon_signed_rank", n, p)
        y = ymax * h
        ax.plot([a, a, b, b], [y, y * 1.01, y * 1.01, y], c="#333", lw=0.8)
        ax.text((a + b) / 2, y * 1.02, f"$p$={fmtp(p)}", ha="center", fontsize=14)
    ax.set_ylim(top=ymax * 1.35)  # keep the p annotations inside the axes frame
    ax.set_xticks(range(3))
    ax.set_xticklabels([LABELS[c] for c in ORDER], fontsize=17)
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.tick_params(axis="y", labelsize=16)
    ax.set_ylabel("Median Jaccard score", fontsize=17)
    ax.set_title("Fine-grained cell-type pair", fontsize=18)
    ax.set_axisbelow(True)
    ax.grid(True, color="#DCDCDC", linewidth=0.8)
    fig.tight_layout()
    savefig(fig, outdir, out_name(name, prefix))


def fig3gh(sheet: pd.DataFrame, anchor: str, tag: str, outdir: Path, source_dir: Path,
           prefix: str = "fig3") -> None:
    found = sheet["anchor"].unique().tolist()
    if found != [anchor]:
        raise SystemExit(f"[{prefix}{tag}] expected anchor {anchor!r} in the sheet, found {found}")
    # partner order = sheet row order (Boundary-median ranking, top-n already applied)
    top = list(dict.fromkeys(sheet["partner"]))
    panel = f"fig3{tag}_{anchor}_partners"
    src = sheet[["anchor", "partner", "sample_id"] + [LABELS[c] for c in ORDER]]
    src.to_csv(source_dir / f"source_data_{prefix}{tag}.csv", index=False)

    fig, axes = plt.subplots(1, len(top), figsize=(1.25 * len(top) + 1.2, 9.0), sharey=True)
    axes = np.atleast_1d(axes)
    for part, ax in zip(top, axes):
        piv = compartment_matrix(sheet[sheet["partner"] == part], "sample_id").reindex(columns=ORDER)
        for i, comp in enumerate(ORDER):
            v = piv[comp].dropna().to_numpy()
            ax.boxplot([v], positions=[i], widths=0.72, patch_artist=True,
                       boxprops={"facecolor": COMP_COLORS[comp], "edgecolor": "#111", "linewidth": 0.9},
                       medianprops={"color": "#111", "linewidth": 1.2},
                       whiskerprops={"color": "#111", "linewidth": 0.9},
                       capprops={"color": "#111", "linewidth": 0.9},
                       flierprops={"marker": "o", "markersize": 2.2, "markerfacecolor": "#111",
                                   "markeredgecolor": "none", "alpha": 0.8})
        for (a, b), y in [((0, 1), 1.05), ((1, 2), 1.18), ((0, 2), 1.31)]:
            p, n = paired_wilcoxon(piv, ORDER[a], ORDER[b])
            record(panel, part, f"{ORDER[a]}-{ORDER[b]}", "wilcoxon_signed_rank", n, p)
            if np.isfinite(p):
                ax.plot([a, a, b, b], [y, y + 0.015, y + 0.015, y], c="#111", lw=0.7)
                ax.text(1.0, y + 0.03, f"$p$={fmtp(p)}", ha="center", va="bottom", fontsize=12.5)
        ax.set_xlim(-0.75, 2.75)
        ax.set_ylim(-0.04, 1.48)
        ax.set_yticks([0.0, 0.5, 1.0])
        ax.set_xticks([1])
        ax.set_xticklabels([part], rotation=45, ha="right", fontsize=22)
        ax.tick_params(axis="y", labelsize=22)
        for s in ax.spines.values():
            s.set_linewidth(0.8)
            s.set_color("#333")
    axes[0].set_ylabel("Jaccard score", fontsize=24)
    fig.suptitle(anchor.replace("_", " "), fontsize=24, x=0.575, y=0.965)
    # fixed margins: fig3g and fig3h share the same canvas size
    fig.subplots_adjust(left=0.195, right=0.995, top=0.92, bottom=0.30, wspace=0.10)
    savefig(fig, outdir, out_name(panel, prefix), tight=False)


def run_panel_set(name: str, panels: list) -> None:
    spec = PANEL_SETS[name]
    outdir, source_dir = cfg.output_dirs(spec["figure"])
    prefix = spec["prefix"]
    STATS.clear()
    if "c" in panels:
        fig3c(cfg.read_source(spec["sheets"]["c"]), outdir, source_dir, prefix)
    for tag, anchor in ANCHORS.items():
        if tag in panels:
            fig3gh(cfg.read_source(spec["sheets"][tag]), anchor, tag, outdir, source_dir, prefix)
    if "g" in panels or "h" in panels:
        compartment_legend(outdir, prefix)
    stats = pd.DataFrame(STATS)
    stats_out = source_dir / f"{prefix}cgh_statistics.csv"
    stats.to_csv(stats_out, index=False)
    print("stats table:", stats_out, f"({len(stats)} tests)")


def main() -> None:
    ap = argparse.ArgumentParser(description="Fig 3c/g/h and Extended Data Fig 7c/g/h from Source Data "
                                             "(paired Wilcoxon signed-rank).",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--panel", choices=["fig3", "extfig7", "all"], default="all",
                    help="fig3 = Fig 3c/g/h (Source_Data_Fig3.xlsx); extfig7 = Extended Data Fig 7c/g/h "
                         "(Source_Data_Extended_Fig7.xlsx, BRCA-downsampled cohort)")
    ap.add_argument("--panels", nargs="+", default=["c", "g", "h"], choices=["c", "g", "h"],
                    help="sub-panels to draw")
    args = ap.parse_args()

    for name in (["fig3", "extfig7"] if args.panel == "all" else [args.panel]):
        run_panel_set(name, args.panels)


if __name__ == "__main__":
    main()
