#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fig. 4f,g and Supp. Figs. 12-13: TCGA Kaplan-Meier curves (HH/HL/LH/LL) and Cox forest plots per LR pair.

Input : Source_Data_Fig4.xlsx, cfg.read_source "Fig4f" (Cox) / "Fig4g" (per-tumour KM), PLAU-PLAUR;
        Source_Data_Supp_Fig.xlsx, cfg.read_source_sections "SuppFig12_top" / "SuppFig12_bottom" /
        "SuppFig13", FN1-PLAUR / SPP1-ITGA5 / SERPINE1-PLAUR.
Output: Plots/survival/<ct_pair>/<lr_pair>/km_4group*.pdf and cox_forest_HHvsLL_viridis_q_per_lr.pdf,
        Tables/lr_ct_cox_summary.csv under Output/Figure4/ (--outdir changes the root).
Only the log-rank P values on the KM plots are recomputed (lifelines 0.30); HR, 95% CI
and BH q are the sheet values.
The per-tumour sheets hold only the cancer types shown, so KM curves are drawn for those. Plotting code is
that of Analysis/Run_Fig4_Statistics/03_Fig4fg_SupFig12-13_survival_cox_km.py, so the PDFs match.
"""
from __future__ import annotations

# "Libraries and paths" ----
import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # Code/ for `common`
from common import config as cfg

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.ticker import FuncFormatter


# "Inputs" ----
# LR pair -> Source Data: a tuple gives (Cox key, KM key) for cfg.read_source; a string names one sheet whose
# 'Cox proportional hazard model' / 'Kaplan-meier' sections are read with cfg.read_source_sections.
SOURCE_PANELS = {
    "PLAU-PLAUR": ("Fig4f", "Fig4g"),      # Fig. 4f (Cox forest), Fig. 4g (KM)
    "FN1-PLAUR": "SuppFig12_top",          # Supp. Fig. 12 top
    "SPP1-ITGA5": "SuppFig12_bottom",      # Supp. Fig. 12 bottom
    "SERPINE1-PLAUR": "SuppFig13",         # Supp. Fig. 13
}
COX_SECTION, KM_SECTION = "Cox proportional hazard model", "Kaplan-meier"
# sheet column -> internal name
COX_COLUMNS = {
    "Cell-type pair": "ct_pair", "Ligand-receptor pair": "lr_pair", "Cancer group": "cancer_group",
    "HR (HH vs LL)": "hr_HH_vs_LL", "95% CI lower": "ci_lower_HH_vs_LL", "95% CI upper": "ci_upper_HH_vs_LL",
    "BH q (per-LR)": "cox_q_HH_vs_LL_per_lr",
}
KM_COLUMNS = {
    "Sample barcode": "sample_barcode", "Cancer group": "cancer_group", "Ligand-receptor pair": "lr_pair",
    "Overall survival time (days)": "OS_time", "Overall survival status (1=death)": "OS_event",
    "Joint group": "joint_group",
}

# "Outputs" ----
# Default root is Output/Figure4/ (override with --outdir):
#   <outdir>/Plots/survival/<ct_pair>/<lr_pair>/km_4group*.pdf, cox_forest_*.pdf
#   <outdir>/Tables/lr_ct_cox_summary.csv
FIG4_PLOTS, FIG4_TABLES = cfg.output_dirs("Fig4")
DEFAULT_OUTDIR = FIG4_PLOTS.parent
SURV_OUT = FIG4_PLOTS / "survival"
SUMMARY_CSV = FIG4_TABLES / "lr_ct_cox_summary.csv"


def set_output_paths(outdir: Path) -> None:
    """Point SURV_OUT / SUMMARY_CSV at <outdir>/Plots/survival and <outdir>/Tables."""
    global SURV_OUT, SUMMARY_CSV
    SURV_OUT = outdir / "Plots" / "survival"
    SUMMARY_CSV = outdir / "Tables" / "lr_ct_cox_summary.csv"


CT_LABEL = "HSPA6+ Macro - tCAF"

LR_PAIRS = ["PLAU-PLAUR", "FN1-PLAUR", "SPP1-ITGA5", "SERPINE1-PLAUR"]

TARGET_CANCER_GROUPS = ["BRCA", "LUCA", "KICA", "COCA", "HNCA", "SKCA", "STCA", "LICA", "PACA"]
P_FLOOR = 1e-12
NEG_LOG_P_VMAX_MIN = 3.0


def configure_plot_style() -> None:
    # Falls back to the default sans-serif where Arial is not installed.
    try:
        font_manager.findfont("Arial", fallback_to_default=False)
        family = "Arial"
    except ValueError:
        family = plt.rcParams["font.sans-serif"][0]
        print(f"[configure_plot_style] Arial not found; using '{family}' instead. "
              "Install Arial to reproduce the submitted figure typography.", flush=True)
    plt.rcParams.update({
        "font.family": family,
        "font.sans-serif": [family],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    if family == "Arial":
        plt.rcParams.update({
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
        })


def sanitize_name(name: str) -> str:
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(name))


def load_source_tables(lr_pair: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Source Data of one LR pair -> (Cox table, per-tumour KM table) with internal column names."""
    spec = SOURCE_PANELS[lr_pair]
    if isinstance(spec, tuple):
        cox, km = cfg.read_source(spec[0]), cfg.read_source(spec[1])
    else:
        sections = cfg.read_source_sections(spec)
        cox, km = sections[COX_SECTION], sections[KM_SECTION]
    for name, df, cols in (("Cox", cox, COX_COLUMNS), ("KM", km, KM_COLUMNS)):
        missing = [c for c in cols if c not in df.columns]
        if missing:
            raise SystemExit(f"[{lr_pair}] {name} Source Data table lacks columns: {missing}")
    cox = cox.rename(columns=COX_COLUMNS)[list(COX_COLUMNS.values())].copy()
    km = km.rename(columns=KM_COLUMNS)[list(KM_COLUMNS.values())].copy()
    if not ((cox["lr_pair"] == lr_pair).all() and (km["lr_pair"] == lr_pair).all()):
        raise SystemExit(f"[{lr_pair}] Source Data table describes another LR pair")
    if not (cox["ct_pair"] == CT_LABEL).all():
        raise SystemExit(f"[{lr_pair}] Source Data Cox table describes another cell-type pair")
    km["OS_event"] = km["OS_event"].astype(int)
    km["OS_time"] = km["OS_time"].astype(float)
    return cox, km


def _fmt_logrank_p(p: float) -> str:
    """Format a log-rank P for the plot: 4 decimals, scientific notation (mathtext) below 1e-4."""
    if p is None or not np.isfinite(p):
        return "n.a."
    if p < 1e-4:
        exp = int(np.floor(np.log10(p)))
        mant = p / 10.0 ** exp
        mant_str = f"{mant:.1f}".rstrip("0").rstrip(".")
        return rf"${mant_str}\times10^{{{exp}}}$"
    return f"{p:.4f}"


def _fmt_q(q: float) -> str:
    """Format an adjusted P for the forest plot: "q=0.0199", scientific notation (mathtext) below 1e-4."""
    if q is None or not np.isfinite(q):
        return "q=n.a."
    if q < 0.0001:
        exp = int(np.floor(np.log10(q)))
        mant = q / 10.0 ** exp
        mant_str = f"{mant:.1f}".rstrip("0").rstrip(".")
        return rf"q={mant_str}$\times10^{{{exp}}}$"
    return f"q={q:.4f}"


def plot_km(gdf: pd.DataFrame, cancer_group: str, lr_pair: str, out_path: Path, km_pdf: PdfPages) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    colors = {"LL": "gray", "LH": "tab:blue", "HL": "tab:orange", "HH": "tab:red"}
    for label in ["LL", "LH", "HL", "HH"]:
        sub = gdf[gdf["joint_group"] == label]
        kmf = KaplanMeierFitter()
        kmf.fit(sub["OS_time"], sub["OS_event"], label=f"{label} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, color=colors[label], ci_show=False)

    hh = gdf[gdf["joint_group"] == "HH"]
    ll = gdf[gdf["joint_group"] == "LL"]
    hl = gdf[gdf["joint_group"] == "HL"]
    lh = gdf[gdf["joint_group"] == "LH"]

    def lrt(a: pd.DataFrame, b: pd.DataFrame):
        return logrank_test(a["OS_time"], b["OS_time"], event_observed_A=a["OS_event"], event_observed_B=b["OS_event"])

    lr_ll = lrt(hh, ll)
    lr_hl = lrt(hh, hl)
    lr_lh = lrt(hh, lh)

    ax.set_title(f"{cancer_group}\n{lr_pair}  ×  {CT_LABEL}", fontsize=18)
    ax.set_xlabel("Overall survival (days)", fontsize=17)
    ax.set_ylabel("Survival probability", fontsize=17)
    _km_lines = [
        "log-rank $P$",
        f"HH vs LL = {_fmt_logrank_p(lr_ll.p_value)}",
        f"HH vs HL = {_fmt_logrank_p(lr_hl.p_value)}",
        f"HH vs LH = {_fmt_logrank_p(lr_lh.p_value)}",
    ]
    # One text per line at a fixed step, so a superscript does not widen its row.
    for _i, _ln in enumerate(reversed(_km_lines)):
        ax.text(0.03, 0.045 + _i * 0.075, _ln, ha="left", va="bottom",
                fontsize=15, transform=ax.transAxes,
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1.2))
    ax.tick_params(axis="both", labelsize=16)
    ax.legend(fontsize=12)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    km_pdf.savefig(fig, dpi=700)
    fig.savefig(out_path, dpi=700, bbox_inches="tight")
    plt.close(fig)


def draw_km_curves(km_tables: dict[str, pd.DataFrame]) -> None:
    """Four-group KM curves per LR pair for the cancer types in the sheet (per-cancer PDFs + multi-page km_4group.pdf)."""
    for lr_pair in LR_PAIRS:
        km = km_tables[lr_pair]
        out_dir = SURV_OUT / sanitize_name(CT_LABEL) / sanitize_name(lr_pair)
        out_dir.mkdir(parents=True, exist_ok=True)
        km_pdf_path = out_dir / "km_4group.pdf"
        drawn = []
        with PdfPages(km_pdf_path) as km_pdf:
            for cancer_group in TARGET_CANCER_GROUPS:
                gdf = km[km["cancer_group"] == cancer_group]
                if gdf.empty:
                    continue        # the sheet holds only the cancer types shown in the figure
                plot_km(gdf, cancer_group, lr_pair, out_dir / f"km_4group_{cancer_group}.pdf", km_pdf)
                drawn.append(f"{cancer_group} (n={len(gdf)})")
        print(f"KM done: {lr_pair}: {', '.join(drawn)}")


def build_summary(cox_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Concatenate the four Cox tables into lr_ct_cox_summary.csv (rows: LR pair order x TARGET_CANCER_GROUPS order)."""
    summary = pd.concat([cox_tables[lr_pair] for lr_pair in LR_PAIRS], ignore_index=True)
    order = pd.DataFrame({
        "_lr": summary["lr_pair"].map({p: i for i, p in enumerate(LR_PAIRS)}),
        "_cg": summary["cancer_group"].map({c: i for i, c in enumerate(TARGET_CANCER_GROUPS)}),
    })
    if order.isna().any().any():
        unknown = summary.loc[order.isna().any(axis=1), ["lr_pair", "cancer_group"]]
        raise SystemExit(f"[Cox] unexpected LR pair / cancer group in Source Data:\n{unknown}")
    summary = summary.loc[order.sort_values(["_lr", "_cg"]).index].reset_index(drop=True)
    SUMMARY_CSV.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(SUMMARY_CSV, index=False)
    return summary


def _p_col_and_label() -> tuple[str, str]:
    return "cox_q_HH_vs_LL_per_lr", r"-log10(adjusted $P$-value)"


def plot_forest(summary: pd.DataFrame) -> None:
    p_col, cbar_label = _p_col_and_label()
    for (ct_pair, lr_pair), pair_df in summary.groupby(["ct_pair", "lr_pair"], sort=False):
        plot_df = pair_df.dropna(subset=["hr_HH_vs_LL", p_col]).sort_values("hr_HH_vs_LL").reset_index(drop=True)
        if plot_df.empty:
            continue
        out_dir = SURV_OUT / sanitize_name(ct_pair) / sanitize_name(lr_pair)
        out_dir.mkdir(parents=True, exist_ok=True)
        n = len(plot_df)
        neg_log_p = -np.log10(plot_df[p_col].clip(lower=P_FLOOR).values)
        norm = Normalize(vmin=0.0, vmax=max(float(neg_log_p.max()), NEG_LOG_P_VMAX_MIN))

        fig, ax = plt.subplots(figsize=(7.5, 0.4 * n + 1.5))
        _max_ci = plot_df["ci_upper_HH_vs_LL"].max()
        _min_ci = plot_df["ci_lower_HH_vs_LL"].min()
        _span_ci = max(_max_ci - _min_ci, 0.01)
        _x_q = _max_ci + _span_ci * 0.03   # q column just right of the widest CI
        for idx, row in plot_df.iterrows():
            color = plt.cm.viridis(norm(-np.log10(max(row[p_col], P_FLOOR))))
            edge = "red" if row[p_col] < 0.05 else "black"
            err_low = row["hr_HH_vs_LL"] - row["ci_lower_HH_vs_LL"]
            err_high = row["ci_upper_HH_vs_LL"] - row["hr_HH_vs_LL"]
            ax.errorbar(row["hr_HH_vs_LL"], idx, xerr=[[err_low], [err_high]],
                        fmt="o", color=color, ecolor="gray", capsize=4,
                        markersize=11, markeredgecolor=edge, markeredgewidth=2 if edge == "red" else 0.6)
            _qv = row[p_col]
            _qs = _fmt_q(_qv)
            ax.text(_x_q, idx, _qs, ha="left", va="center", fontsize=15, color="black")
        ax.axvline(1.0, color="black", linestyle="--", linewidth=1)
        min_ci = plot_df["ci_lower_HH_vs_LL"].min()
        max_ci = plot_df["ci_upper_HH_vs_LL"].max()
        span = max(max_ci - min_ci, 0.01)
        ax.set_xlim(max(0.01, min_ci - max(0.05, span * 0.15)), _x_q + _span_ci * 0.28)
        ax.set_ylim(-0.7, n - 1 + 0.75)
        ax.set_yticks(np.arange(n))
        ax.set_yticklabels(plot_df["cancer_group"], fontsize=18)
        ax.set_xlabel("Hazard ratio (HH vs LL)", fontsize=18)
        ax.tick_params(axis="x", labelsize=18)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.2f}"))
        ax.grid(axis="x", alpha=0.3)
        out_path = out_dir / "cox_forest_HHvsLL_viridis_q_per_lr.pdf"
        title_size = 16
        cbar_size = 14

        ax.set_title(f"{lr_pair} × {ct_pair}", fontsize=title_size)
        sm = ScalarMappable(norm=norm, cmap=plt.cm.viridis)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.005, fraction=0.045, shrink=0.6, anchor=(0.0, 0.0))
        cbar.set_label(cbar_label, fontsize=cbar_size)
        cbar.ax.tick_params(labelsize=cbar_size)
        fig.tight_layout()
        fig.savefig(out_path, dpi=700, bbox_inches="tight")
        plt.close(fig)
    print("Forest plots done.")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="TCGA survival figures (KM + Cox forest) of prioritized LR pairs x HSPA6+ Macro-tCAF "
                    "from the paper's Source Data (Fig. 4f,g; Supp. Figs. 12-13).")
    ap.add_argument("--outdir", default=str(DEFAULT_OUTDIR),
                    help="output root: PDFs under <outdir>/Plots/survival/, lr_ct_cox_summary.csv under "
                         "<outdir>/Tables/ (default: %(default)s)")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    set_output_paths(Path(args.outdir).resolve())
    configure_plot_style()
    SURV_OUT.mkdir(parents=True, exist_ok=True)
    print(f"Output: {SURV_OUT}")
    print(f"Cox summary: {SUMMARY_CSV}")
    cox_tables, km_tables = {}, {}
    for lr_pair in LR_PAIRS:
        cox_tables[lr_pair], km_tables[lr_pair] = load_source_tables(lr_pair)
        print(f"Source Data {lr_pair}: Cox rows={len(cox_tables[lr_pair])} | KM rows={len(km_tables[lr_pair])} "
              f"({', '.join(km_tables[lr_pair]['cancer_group'].unique())})")
    draw_km_curves(km_tables)
    summary = build_summary(cox_tables)
    plot_forest(summary)
    print("Done.")


if __name__ == "__main__":
    main()
