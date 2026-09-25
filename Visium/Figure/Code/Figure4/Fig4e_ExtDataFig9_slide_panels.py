#!/usr/bin/env python3
"""Fig. 4e and Ext. Data Fig. 9: per-slide 4 x 2 panels on H&E, compartment, pEMT score
and six co-localization pairs.

Input : Visium AnnData with H&E images (env PCASSO_VISIUM_H5AD or --h5ad-path; not bundled) and the bundled
        per-spot pEMT table Input/Fig4e/pemt_scores.csv.gz (--pemt-csv, --pemt-col).
Output: <outdir>/<cancer_type>/<sample_id>__pemt_overlap_panels.pdf, default outdir
        Output/Figure4/Plots/Fig4e_ExtDataFig9_slide_panels/. Fig. 4e = BRCA/GSM6433585_092A; the other nine
        slides are tiled into Ext. Data Fig. 9 by ExtDataFig9_grid_3x3.py.
Co-localization panels: each feature is arcsinh-transformed (--transform), clipped at --clip-percentile and
scaled to [0, 1]; red = left only, blue = right only, green = overlap (min of the two, raised to
--overlap-gamma), brightness = max of the two (raised to --value-gamma), with a 2-D colour key per panel.
Pairs: tCAF - HSPA6+ Macro, PLAU-PLAUR, FN1-PLAUR, PLAU-ITGB2, SPP1-ITGA5 and a sixth chosen by --pair-set
(serpine1_plaur default, or apoe_trem2).
Main options: --target CANCER_TYPE:SAMPLE_ID[:PEMT_VMAX[:SPOT_SIZE[:VALUE_GAMMA]]] (repeatable; default = the
        ten manuscript slides in two passes), --cancer-label PDAC=PACA (title only; folders keep the obs
        value), --title-format, --row-gap/--panel-height/--col-gap (fixed-row layout, inches), --font-family
        (Arial, embedded as TrueType); see --help.
"""
from __future__ import annotations

# "Libraries and paths" ----
import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # Code/ for `common`
from common import config as cfg

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("MKL_THREADING_LAYER", "GNU")
os.environ.setdefault("KMP_AFFINITY", "disabled")
os.environ.setdefault("KMP_INIT_AT_FORK", "FALSE")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


COMPARTMENT_COLORS = {
    "Malignant": "tab:red",
    "Boundary": "tab:green",
    "Normal": "tab:blue",
}

# "Inputs" ----
# Visium AnnData with H&E images (not bundled; PCASSO_VISIUM_H5AD or --h5ad-path) and the
# bundled per-spot pEMT score table.
DEFAULT_H5AD = cfg.VISIUM_H5AD
DEFAULT_PEMT_CSV = cfg.FIG4E_INPUT / "pemt_scores.csv.gz"

# "Outputs" ----
# One PDF per slide under <outdir>/<cancer_type>/. Fig. 4e = BRCA/GSM6433585_092A;
# the other nine slides are tiled into Extended Data Fig. 9 by ExtDataFig9_grid_3x3.py.
DEFAULT_OUTDIR = cfg.output_dirs("Fig4")[0] / "Fig4e_ExtDataFig9_slide_panels"

# Slides of the manuscript, drawn in two passes when no --target is given:
# pass 1 = the nine slides of Extended Data Fig. 9, pass 2 = the BRCA slide of Fig. 4e.
MANUSCRIPT_PASSES = [
    {"title_format": "{label}", "title_size": 28.0, "col_gap": 2.0,
     "targets": [
         "COCA:scCRLM_Atlas_ST-colon1:0.85",
         "KICA:GSM5924042_frozen_a_1:0.75:27",
         "LUCA:CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_220714:0.7",
         "OVCA:GSM6506110_SP1:0.8",
         "PDAC:GSM7498812_SS1923404:0.9:20:1.4",
         "PRCA:GSM7211257_EHU-W3:0.6",
         "STCA:GSE251950_20_00331_LI_SING:0.9",
         "THCA:S20_63981:0.8",
         "UCEC:01_034_C3d1:0.8",
     ]},
    {"title_format": "{label} | {sample}", "title_size": 22.0, "col_gap": 2.6,
     "targets": ["BRCA:GSM6433585_092A:0.8:28"]},
]

# Co-localization pairs, drawn in this order (left = red, right = blue, overlap = green).
COMMON_OVERLAP_SPECS = [
    {"left": "tCAF", "right": "HSPA6+_Macro",
     "left_label": "tCAF", "right_label": "HSPA6+ Macro",
     "title": "tCAF - HSPA6+ Macro"},
    {"left": "PLAU", "right": "PLAUR",
     "left_label": "PLAU", "right_label": "PLAUR", "title": "PLAU-PLAUR"},
    {"left": "FN1", "right": "PLAUR",
     "left_label": "FN1", "right_label": "PLAUR", "title": "FN1-PLAUR"},
    {"left": "PLAU", "right": "ITGB2",
     "left_label": "PLAU", "right_label": "ITGB2", "title": "PLAU-ITGB2"},
    {"left": "SPP1", "right": "ITGA5",
     "left_label": "SPP1", "right_label": "ITGA5", "title": "SPP1-ITGA5"},
]

# Sixth pair, selected with --pair-set.
SIXTH_PAIR_BY_SET = {
    "apoe_trem2": {"left": "APOE", "right": "TREM2",
                   "left_label": "APOE", "right_label": "TREM2", "title": "APOE-TREM2"},
    "serpine1_plaur": {"left": "SERPINE1", "right": "PLAUR",
                       "left_label": "SERPINE1", "right_label": "PLAUR", "title": "SERPINE1-PLAUR"},
}


def get_overlap_specs(pair_set: str) -> List[Dict[str, str]]:
    if pair_set not in SIXTH_PAIR_BY_SET:
        raise ValueError(f"Unknown --pair-set `{pair_set}`. Choose from {sorted(SIXTH_PAIR_BY_SET)}")
    return COMMON_OVERLAP_SPECS + [SIXTH_PAIR_BY_SET[pair_set]]

COLOR_LEFT_RED = np.array([0.90, 0.12, 0.14], dtype=float)
COLOR_OVERLAP_GREEN = np.array([0.10, 0.75, 0.20], dtype=float)
COLOR_RIGHT_BLUE = np.array([0.10, 0.35, 0.95], dtype=float)


# ---------------------------------------------------------------------------
# Data selection helpers
# ---------------------------------------------------------------------------
def safe_name(s: str) -> str:
    return re.sub(r"[^0-9A-Za-z._+-]+", "_", str(s))


def resolve_slide_col(adata, slide_col: str) -> str:
    if slide_col != "auto":
        if slide_col not in adata.obs.columns:
            raise KeyError(f"`{slide_col}` not in adata.obs")
        return slide_col

    for cand in ("sample_id", "slide_id", "sample"):
        if cand in adata.obs.columns:
            return cand
    raise KeyError("Could not auto-detect slide column. Use --slide-col explicitly.")


def normalize_cancer_type_name(name: str) -> str:
    text = str(name).strip()
    aliases = {"PACA": "PDAC", "UECA": "UCEC", "LIHC": "LICA"}
    return aliases.get(text, text)


def resolve_cancer_type(adata, cancer_col: str, requested: str) -> str:
    if cancer_col not in adata.obs.columns:
        raise KeyError(f"`{cancer_col}` not in adata.obs")
    available = sorted(adata.obs[cancer_col].astype(str).unique().tolist())
    available_norm = {normalize_cancer_type_name(x): x for x in available}
    norm = normalize_cancer_type_name(requested)
    if norm in available_norm:
        return available_norm[norm]
    if requested in available:
        return requested
    raise ValueError(f"cancer_type `{requested}` not found. Available: {available}")


def check_sample_in_cancer_type(adata, cancer_col: str, cancer_type: str, slide_col: str, sample_id: str) -> None:
    mask = adata.obs[cancer_col].astype(str) == str(cancer_type)
    available = set(adata.obs.loc[mask, slide_col].astype(str).unique().tolist())
    if str(sample_id) not in available:
        raise ValueError(f"sample_id `{sample_id}` is not a slide of cancer_type={cancer_type}")


def select_slide(adata, slide_value: str, slide_col: str):
    sub = adata[adata.obs[slide_col].astype(str) == str(slide_value), :].copy()
    if sub.n_obs == 0:
        raise ValueError(f"No spots found for {slide_col}={slide_value}")

    keys = list(sub.uns.get("spatial", {}).keys())
    if not keys:
        raise ValueError("No `uns['spatial']` keys found in selected slide.")

    matched = [k for k in keys if str(slide_value) == str(k) or str(slide_value) in str(k)]
    if len(matched) == 1:
        keep = matched[0]
    elif len(keys) == 1:
        keep = keys[0]
    else:
        keep = matched[0] if matched else keys[0]

    sub.uns["spatial"] = {keep: sub.uns["spatial"][keep]}
    return sub


def attach_pemt_score(adata, pemt_csv: str, pemt_col: str) -> None:
    pemt = pd.read_csv(pemt_csv, index_col=0)
    if pemt_col not in pemt.columns:
        raise KeyError(f"`{pemt_col}` not in pEMT CSV columns: {list(pemt.columns)}")
    adata.obs[pemt_col] = pd.to_numeric(pemt[pemt_col], errors="coerce").reindex(adata.obs_names)


def prepare_compartment_column(
    adata,
    source_col: str = "LocationNew",
    target_col: str = "LocationNew_plot",
) -> None:
    if source_col not in adata.obs.columns:
        raise KeyError(f"`{source_col}` not in adata.obs")
    mapped = adata.obs[source_col].astype(str).replace(
        {
            "Mal": "Malignant",
            "Bdy": "Boundary",
            "Normal": "Normal",
        }
    )
    categories = ["Malignant", "Boundary", "Normal"]
    adata.obs[target_col] = pd.Categorical(mapped, categories=categories, ordered=True)


# ---------------------------------------------------------------------------
# Feature scaling and colour mixing
# ---------------------------------------------------------------------------
def _to_1d(x) -> np.ndarray:
    if sp.issparse(x):
        x = x.toarray()
    return np.asarray(x).ravel()


def _normalize_vector(v: np.ndarray, transform: str, clip_percentile: float) -> np.ndarray:
    """Transform, clip at the given percentile and scale to [0, 1]."""
    out = np.nan_to_num(_to_1d(v).astype(float), nan=0.0, posinf=0.0, neginf=0.0)
    if transform == "arcsinh":
        out = np.arcsinh(out)
    elif transform == "log1p":
        out = np.log1p(out)
    elif transform == "none":
        pass
    else:
        raise ValueError(f"Unsupported transform: {transform}")

    hi = np.nanpercentile(out, clip_percentile)
    if hi <= 0:
        return np.zeros_like(out)
    return np.clip(out, 0, hi) / hi


def _get_feature_vector(adata, name: str, layer: str | None, use_raw: bool) -> np.ndarray:
    if name in adata.obs.columns:
        return adata.obs[name].to_numpy()

    if name in adata.var_names:
        if layer is None:
            return adata[:, name].X
        idx = int(adata.var_names.get_loc(name))
        return adata.layers[layer][:, idx]

    if use_raw and hasattr(adata, "raw") and adata.raw is not None and name in adata.raw.var_names:
        return adata.raw[:, name].X

    raise KeyError(f"Feature `{name}` not found in adata.obs/var_names/raw.var_names.")


def _mix_overlap_green_rgb(
    left_v: np.ndarray,
    right_v: np.ndarray,
    overlap_gamma: float,
    value_gamma: float,
) -> np.ndarray:
    """Mix two [0, 1] vectors into RGB: red = left only, blue = right only, green = overlap (min, ^overlap_gamma);
    brightness = max of the two, ^value_gamma (exponents below 1 widen the green region / brighten dim spots)."""
    left_v = np.clip(np.asarray(left_v, dtype=float).ravel(), 0.0, 1.0)
    right_v = np.clip(np.asarray(right_v, dtype=float).ravel(), 0.0, 1.0)

    overlap = np.minimum(left_v, right_v)
    overlap_w = np.power(overlap, overlap_gamma) if overlap_gamma > 0 and overlap_gamma != 1.0 else overlap
    left_only = np.clip(left_v - overlap, 0.0, 1.0)
    right_only = np.clip(right_v - overlap, 0.0, 1.0)

    w_red = left_only
    w_green = overlap_w
    w_blue = right_only
    w_sum = w_red + w_green + w_blue

    rgb = np.zeros((left_v.size, 3), dtype=float)
    valid = w_sum > 1e-12
    rgb[valid] = (
        w_red[valid, None] * COLOR_LEFT_RED[None, :]
        + w_green[valid, None] * COLOR_OVERLAP_GREEN[None, :]
        + w_blue[valid, None] * COLOR_RIGHT_BLUE[None, :]
    ) / w_sum[valid, None]

    value = np.maximum(left_v, right_v)
    if value_gamma > 0 and value_gamma != 1.0:
        value = np.power(value, value_gamma)
    rgb = rgb * value[:, None]
    return np.clip(rgb, 0.0, 1.0)


def _pair_to_rgb(
    adata,
    left_name: str,
    right_name: str,
    transform: str,
    clip_percentile: float,
    layer: str | None,
    use_raw: bool,
    overlap_gamma: float,
    value_gamma: float,
) -> np.ndarray:
    left_v = _normalize_vector(
        _get_feature_vector(adata, left_name, layer=layer, use_raw=use_raw),
        transform=transform,
        clip_percentile=clip_percentile,
    )
    right_v = _normalize_vector(
        _get_feature_vector(adata, right_name, layer=layer, use_raw=use_raw),
        transform=transform,
        clip_percentile=clip_percentile,
    )
    return _mix_overlap_green_rgb(
        left_v=left_v,
        right_v=right_v,
        overlap_gamma=overlap_gamma,
        value_gamma=value_gamma,
    )


# ---------------------------------------------------------------------------
# Spatial image helpers
# ---------------------------------------------------------------------------
def _get_spatial_image_and_xy(slide, img_key: str):
    lib = list(slide.uns["spatial"].keys())[0]
    img = slide.uns["spatial"][lib]["images"][img_key]
    sf_key = f"tissue_{img_key}_scalef"
    sf = slide.uns["spatial"][lib]["scalefactors"][sf_key]
    xy = slide.obsm["spatial"] * sf
    return img, xy


def compute_spot_crop_box(slide, img_key: str, pad_frac: float = 0.05):
    """Bounding box of the spots plus padding, clipped to the image: (x0, x1, y0, y1) or None."""
    img, xy = _get_spatial_image_and_xy(slide, img_key)
    x = xy[:, 0]
    y = xy[:, 1]
    finite = np.isfinite(x) & np.isfinite(y)
    if not finite.any():
        return None
    x = x[finite]
    y = y[finite]
    x0, x1 = float(x.min()), float(x.max())
    y0, y1 = float(y.min()), float(y.max())
    pad_x = (x1 - x0) * pad_frac
    pad_y = (y1 - y0) * pad_frac
    x0 -= pad_x
    x1 += pad_x
    y0 -= pad_y
    y1 += pad_y
    h, w = int(img.shape[0]), int(img.shape[1])
    return max(x0, 0.0), min(x1, float(w)), max(y0, 0.0), min(y1, float(h))


def apply_spot_crop(ax, slide, img_key: str, pad_frac: float = 0.05) -> None:
    """Restrict the axes view to the spot bounding box plus padding.
    Insets (colourbar, legend) are placed in axes coordinates, so the crop does not move them."""
    box = compute_spot_crop_box(slide, img_key, pad_frac)
    if box is None:
        return
    x0, x1, y0, y1 = box
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)  # image origin is upper-left, so invert y


# ---------------------------------------------------------------------------
# Panel drawing
# ---------------------------------------------------------------------------
def draw_overlap_colorbar_right(
    ax,
    left_label: str,
    right_label: str,
    overlap_gamma: float,
    value_gamma: float,
    tick_size: float,
) -> None:
    """2-D colour legend (left feature on x, right feature on y) to the right of the panel."""
    inset = ax.inset_axes([1.20, 0.12, 0.18, 0.72], transform=ax.transAxes)
    grid_size = 100
    x = np.linspace(0.0, 1.0, grid_size, dtype=float)
    y = np.linspace(0.0, 1.0, grid_size, dtype=float)
    left_grid, right_grid = np.meshgrid(x, y)
    rgb = _mix_overlap_green_rgb(
        left_v=left_grid.ravel(),
        right_v=right_grid.ravel(),
        overlap_gamma=overlap_gamma,
        value_gamma=value_gamma,
    ).reshape(grid_size, grid_size, 3)

    inset.imshow(rgb, origin="lower", extent=[0, 1, 0, 1], interpolation="nearest")
    inset.set_xlabel(left_label, fontsize=tick_size, labelpad=1)
    inset.set_ylabel(right_label, fontsize=tick_size, labelpad=1)
    inset.set_xticks([0.0, 1.0])
    inset.set_yticks([0.0, 0.5, 1.0])
    inset.tick_params(axis="both", labelsize=tick_size, length=1.5, pad=1)
    for spine in inset.spines.values():
        spine.set_linewidth(0.6)


def plot_compartment_panel(
    ax,
    slide,
    img_key: str,
    title: str,
    spot_size: float,
    alpha_img: float,
    alpha_spot: float,
    title_size: float,
    legend_fontsize: float,
) -> None:
    img, xy = _get_spatial_image_and_xy(slide, img_key)
    ax.imshow(img, alpha=alpha_img)
    cats = slide.obs["LocationNew_plot"].astype(str).to_numpy()
    categories = ["Malignant", "Boundary", "Normal"]
    for cat in categories:
        mask = cats == cat
        if mask.any():
            ax.scatter(
                xy[mask, 0],
                xy[mask, 1],
                s=spot_size,
                c=COMPARTMENT_COLORS[cat],
                edgecolors="none",
                alpha=alpha_spot,
                label=cat,
            )
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=COMPARTMENT_COLORS[cat], markersize=6, label=cat)
        for cat in categories
    ]
    ax.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.95, 0.33),
        bbox_transform=ax.transAxes,
        frameon=False,
        fontsize=legend_fontsize,
        handletextpad=0.4,
        borderpad=0.2,
        markerscale=1.5,
    )
    ax.set_title(title, fontsize=title_size)
    ax.set_axis_off()


def plot_pemt_panel(
    ax,
    slide,
    pemt_col: str,
    img_key: str,
    title: str,
    spot_size: float,
    alpha_img: float,
    alpha_spot: float,
    cmap: str,
    pemt_vmax: float,
    title_size: float,
    colorbar_tick_size: float,
    pemt_label: str = "Score",
    pemt_label_size: float | None = None,
) -> None:
    img, xy = _get_spatial_image_and_xy(slide, img_key)
    ax.imshow(img, alpha=alpha_img)
    values = pd.to_numeric(slide.obs[pemt_col], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(values)
    if (~valid).any():
        ax.scatter(
            xy[~valid, 0],
            xy[~valid, 1],
            s=spot_size,
            c="lightgray",
            edgecolors="none",
            alpha=alpha_spot,
        )
    if valid.any():
        sc_obj = ax.scatter(
            xy[valid, 0],
            xy[valid, 1],
            s=spot_size,
            c=values[valid],
            cmap=cmap,
            vmin=0.0,
            vmax=pemt_vmax,
            edgecolors="none",
            alpha=alpha_spot,
        )
        cbar_ax = ax.inset_axes([1.0, 0.005, 0.035, 0.5], transform=ax.transAxes)
        cb = plt.colorbar(sc_obj, cax=cbar_ax)
        if pemt_label:
            cb.set_label(
                pemt_label,
                fontsize=(pemt_label_size if pemt_label_size is not None else colorbar_tick_size),
                labelpad=4,
            )
        cbar_ax.tick_params(axis="both", labelsize=colorbar_tick_size, length=1.5, pad=1)
        for spine in cbar_ax.spines.values():
            spine.set_linewidth(0.6)
    ax.set_title(title, fontsize=title_size)
    ax.set_axis_off()


def plot_overlap_panel(
    ax,
    slide,
    left_name: str,
    right_name: str,
    left_label: str,
    right_label: str,
    title: str,
    img_key: str,
    spot_size: float,
    alpha_img: float,
    transform: str,
    clip_percentile: float,
    layer: str | None,
    use_raw: bool,
    overlap_gamma: float,
    value_gamma: float,
    title_size: float,
    legend_tick_size: float,
) -> None:
    img, xy = _get_spatial_image_and_xy(slide, img_key)
    rgb = _pair_to_rgb(
        slide,
        left_name=left_name,
        right_name=right_name,
        transform=transform,
        clip_percentile=clip_percentile,
        layer=layer,
        use_raw=use_raw,
        overlap_gamma=overlap_gamma,
        value_gamma=value_gamma,
    )
    ax.imshow(img, alpha=alpha_img)
    ax.scatter(xy[:, 0], xy[:, 1], s=spot_size, c=rgb, edgecolors="none")
    ax.set_title(title, fontsize=title_size)
    ax.set_axis_off()
    draw_overlap_colorbar_right(
        ax=ax,
        left_label=left_label,
        right_label=right_label,
        overlap_gamma=overlap_gamma,
        value_gamma=value_gamma,
        tick_size=legend_tick_size,
    )


def plot_sample_panel(
    slide,
    sample_id: str,
    cancer_type: str,
    out_pdf: Path,
    overlap_specs: Sequence[Dict[str, str]],
    pemt_col: str,
    pemt_vmax: float,
    img_key: str,
    spot_size: float,
    alpha_img: float,
    alpha_spot: float,
    cmap: str,
    dpi: int,
    panel_title_size: float,
    colorbar_tick_size: float,
    legend_fontsize: float,
    transform: str,
    clip_percentile: float,
    layer: str | None,
    use_raw: bool,
    overlap_gamma: float,
    value_gamma: float,
    pemt_label: str = "Score",
    pemt_label_size: float | None = None,
    crop_to_spots: bool = True,
    crop_pad_frac: float = 0.05,
    title_gap: float = 0.15,
    title_size: float = 19,
    title_label: str | None = None,
    title_format: str = "{label} | {sample}",
    row_gap: float | None = None,
    panel_height: float | None = None,
    col_gap: float = 2.0,
) -> None:
    fig_w = 16.0
    if row_gap is None:
        # Default layout: fixed 16 x 29 in figure, spacing settled by tight_layout.
        fig, axes = plt.subplots(4, 2, figsize=(fig_w, 29),
                                 gridspec_kw={"wspace": 0.3, "hspace": 0.30})
        fixed_layout = False
    else:
        # Fixed-row layout: every grid cell has the aspect ratio of the cropped image, so no blank space
        # is left above/below the panels; rows sit `row_gap` inches apart (the gap holds the panel titles).
        box = compute_spot_crop_box(slide, img_key, crop_pad_frac) if crop_to_spots else None
        if box is None:
            img, _ = _get_spatial_image_and_xy(slide, img_key)
            aspect = img.shape[0] / img.shape[1]
        else:
            x0, x1, y0, y1 = box
            aspect = (y1 - y0) / max(x1 - x0, 1e-9)
        # The column gap is in inches: the legends right of each panel have a fixed size in points,
        # and a gap proportional to panel width would swallow them on narrow (tall-tissue) slides.
        left, right = 0.02, 0.98
        if panel_height is None:
            panel_w = ((right - left) * fig_w - col_gap) / 2
            panel_h = panel_w * aspect
        else:
            panel_h = float(panel_height)
            panel_w = panel_h / max(aspect, 1e-9)
            fig_w = (2 * panel_w + col_gap) / (right - left)
        top_margin = 1.0      # figure title + first-row panel titles
        bottom_margin = 0.25
        fig_h = top_margin + 4 * panel_h + 3 * row_gap + bottom_margin
        fig, axes = plt.subplots(4, 2, figsize=(fig_w, fig_h))
        fig.subplots_adjust(left=left, right=right,
                            top=1 - top_margin / fig_h, bottom=bottom_margin / fig_h,
                            wspace=col_gap / panel_w, hspace=row_gap / panel_h)
        fixed_layout = True
    axes_flat = axes.ravel()

    plot_compartment_panel(
        ax=axes_flat[0],
        slide=slide,
        img_key=img_key,
        title="Compartment",
        spot_size=spot_size,
        alpha_img=alpha_img,
        alpha_spot=alpha_spot,
        title_size=panel_title_size,
        legend_fontsize=legend_fontsize,
    )

    plot_pemt_panel(
        ax=axes_flat[1],
        slide=slide,
        pemt_col=pemt_col,
        img_key=img_key,
        title="pEMT Signature",
        spot_size=spot_size,
        alpha_img=alpha_img,
        alpha_spot=alpha_spot,
        cmap=cmap,
        pemt_vmax=pemt_vmax,
        title_size=panel_title_size,
        colorbar_tick_size=colorbar_tick_size,
        pemt_label=pemt_label,
        pemt_label_size=pemt_label_size,
    )

    for ax, spec in zip(axes_flat[2:], overlap_specs):
        plot_overlap_panel(
            ax=ax,
            slide=slide,
            left_name=spec["left"],
            right_name=spec["right"],
            left_label=spec["left_label"],
            right_label=spec["right_label"],
            title=spec["title"],
            img_key=img_key,
            spot_size=spot_size,
            alpha_img=alpha_img,
            transform=transform,
            clip_percentile=clip_percentile,
            layer=layer,
            use_raw=use_raw,
            overlap_gamma=overlap_gamma,
            value_gamma=value_gamma,
            title_size=panel_title_size,
            legend_tick_size=max(7.0, colorbar_tick_size - 1),
        )

    n_used = 2 + len(overlap_specs)
    for ax in axes_flat[n_used:]:
        ax.set_axis_off()

    if crop_to_spots:
        for ax in axes_flat[:n_used]:
            apply_spot_crop(ax, slide, img_key, pad_frac=crop_pad_frac)

    if not fixed_layout:
        fig.tight_layout(rect=(0, 0, 1, 0.985))

    # Figure title just above the top row; axes can be shorter than their grid cells (images keep
    # their aspect ratio), so the top edge is measured after layout rather than assumed.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_fig = fig.transFigure.inverted()
    top = max(ax.get_tightbbox(renderer).transformed(to_fig).y1 for ax in axes_flat[:2])
    fig.suptitle(
        title_format.format(label=title_label or cancer_type, cancer=cancer_type, sample=sample_id),
        fontsize=title_size,
        y=top + title_gap / fig.get_figheight(),
        va="bottom",
    )
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, dpi=dpi, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_target(text: str) -> Dict[str, object]:
    """Parse CANCER_TYPE:SAMPLE_ID[:PEMT_VMAX[:SPOT_SIZE[:VALUE_GAMMA]]]."""
    parts = [p.strip() for p in str(text).split(":")]
    if len(parts) not in (2, 3, 4, 5) or not parts[0] or not parts[1]:
        raise argparse.ArgumentTypeError(
            f"--target must be CANCER_TYPE:SAMPLE_ID[:PEMT_VMAX[:SPOT_SIZE[:VALUE_GAMMA]]], got `{text}`"
        )

    def _num(value: str, label: str) -> float | None:
        if not value:
            return None
        try:
            return float(value)
        except ValueError as e:
            raise argparse.ArgumentTypeError(f"{label} must be numeric in `{text}`") from e

    vmax = _num(parts[2], "PEMT_VMAX") if len(parts) >= 3 else None
    spot = _num(parts[3], "SPOT_SIZE") if len(parts) >= 4 else None
    vgamma = _num(parts[4], "VALUE_GAMMA") if len(parts) == 5 else None
    return {"cancer_type": parts[0], "sample_id": parts[1], "pemt_vmax": vmax,
            "spot_size": spot, "value_gamma": vgamma}


def parse_label(text: str) -> tuple[str, str]:
    """Parse CANCER_TYPE=LABEL."""
    if "=" not in text:
        raise argparse.ArgumentTypeError(f"--cancer-label must be CANCER_TYPE=LABEL, got `{text}`")
    key, label = (p.strip() for p in text.split("=", 1))
    if not key or not label:
        raise argparse.ArgumentTypeError(f"--cancer-label must be CANCER_TYPE=LABEL, got `{text}`")
    return key, label


def configure_fonts(font_family: str) -> None:
    """Use one sans-serif family for all text and embed fonts as TrueType in PDF/PS."""
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [font_family, "Liberation Sans", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fig. 4e: compartment / pEMT / co-localization spatial panels per Visium slide.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--h5ad-path", default=str(DEFAULT_H5AD),
                        help="Concatenated Visium AnnData (.h5ad); not bundled, see PCASSO_VISIUM_H5AD.")
    parser.add_argument("--pemt-csv", default=str(DEFAULT_PEMT_CSV),
                        help="Per-spot pEMT score CSV (index = spot barcode).")
    parser.add_argument("--pemt-col", default="pEMT_Itaiyanai", help="pEMT score column in --pemt-csv.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory.")
    parser.add_argument("--target", action="append", type=parse_target, default=None,
                        help="CANCER_TYPE:SAMPLE_ID[:PEMT_VMAX]; repeat once per slide. "
                             "Default: the ten manuscript slides, drawn in two passes "
                             "(the Extended Data Fig. 9 set, then the Fig. 4e BRCA slide).")
    parser.add_argument("--pemt-vmax", type=float, default=0.9,
                        help="Default upper colour limit of the pEMT panel (used when a target omits it).")
    parser.add_argument("--cancer-label", action="append", type=parse_label, default=[],
                        help="CANCER_TYPE=LABEL; title label for that cancer type (e.g. PDAC=PACA). Repeatable.")
    parser.add_argument("--pair-set", choices=sorted(SIXTH_PAIR_BY_SET), default="serpine1_plaur",
                        help="Sixth co-localization pair: apoe_trem2 (APOE vs TREM2) or "
                             "serpine1_plaur (SERPINE1 vs PLAUR).")
    parser.add_argument("--cancer-col", default="cancer_type")
    parser.add_argument("--slide-col", default="auto")
    parser.add_argument("--img-key", default="hires")
    parser.add_argument("--spot-size", type=float, default=20.0,
                        help="Default scatter marker size (used when a target omits SPOT_SIZE).")
    parser.add_argument("--alpha-img", type=float, default=1.0)
    parser.add_argument("--alpha-spot", type=float, default=0.9)
    parser.add_argument("--cmap", default="viridis", help="Colormap of the pEMT panel.")
    parser.add_argument("--dpi", type=int, default=700)
    parser.add_argument("--panel-title-size", type=float, default=24)
    parser.add_argument("--colorbar-tick-size", type=float, default=20)
    parser.add_argument("--legend-fontsize", type=float, default=20)
    parser.add_argument("--transform", choices=["arcsinh", "log1p", "none"], default="arcsinh",
                        help="Transform applied to each feature before scaling.")
    parser.add_argument("--clip-percentile", type=float, default=95.0,
                        help="Percentile used as the upper limit when scaling each feature to [0, 1].")
    parser.add_argument("--overlap-gamma", type=float, default=1.3,
                        help="Exponent on the overlap (green) weight; <1 widens, >1 narrows the green region.")
    parser.add_argument("--value-gamma", type=float, default=1.0,
                        help="Exponent on spot brightness; <1 brightens low/intermediate spots.")
    parser.add_argument("--pemt-label", default="Score", help="pEMT colourbar label (empty string to hide).")
    parser.add_argument("--pemt-label-size", type=float, default=None,
                        help="pEMT colourbar label font size (defaults to --colorbar-tick-size).")
    parser.add_argument("--layer", default=None, help="AnnData layer to read genes from (default: X).")
    parser.add_argument("--use-raw", action="store_true", help="Fall back to adata.raw for missing genes.")
    parser.add_argument("--crop-to-spots", dest="crop_to_spots", action="store_true", default=True,
                        help="Crop each panel to the spot-covered area (default).")
    parser.add_argument("--no-crop", dest="crop_to_spots", action="store_false",
                        help="Show the full H&E image without cropping.")
    parser.add_argument("--crop-pad-frac", type=float, default=0.05,
                        help="Padding around the spot bounding box, as a fraction of its size.")
    parser.add_argument("--title-gap", type=float, default=0.15,
                        help="Gap between the figure title (cancer type | sample) and the top row, in inches.")
    parser.add_argument("--title-size", type=float, default=22, help="Figure title font size.")
    parser.add_argument("--title-format", default="{label} | {sample}",
                        help="Figure title template; fields: {label} (cancer label, see --cancer-label), "
                             "{cancer} (obs cancer type), {sample} (sample id). E.g. \"{label}\" for the type only.")
    parser.add_argument("--row-gap", type=float, default=0.5,
                        help="Fixed-row layout: vertical gap between panel rows in inches (holds the panel titles).")
    parser.add_argument("--panel-height", type=float, default=4.6,
                        help="Fixed-row layout: panel height in inches (figure width follows). Default fills a 16 in wide figure.")
    parser.add_argument("--col-gap", type=float, default=2.0,
                        help="Fixed-row layout: gap between the two panel columns in inches (holds the legends).")
    parser.add_argument("--font-family", default="Arial",
                        help="Font family for all text in the figure (falls back to Liberation Sans / DejaVu Sans).")
    args = parser.parse_args()
    cfg.require_external(Path(args.h5ad_path), "Visium AnnData (--h5ad-path)")
    overlap_specs = get_overlap_specs(args.pair_set)
    configure_fonts(args.font_family)
    cancer_labels = {"PDAC": "PACA", "UCEC": "UECA", **dict(args.cancer_label)}

    adata = sc.read_h5ad(args.h5ad_path)
    slide_col = resolve_slide_col(adata, args.slide_col)

    location_col = next((c for c in ("LocationNew", "Location") if c in adata.obs.columns), None)
    if location_col is None:
        raise KeyError("`LocationNew` (or `Location`) is required in adata.obs")

    attach_pemt_score(adata=adata, pemt_csv=args.pemt_csv, pemt_col=args.pemt_col)
    prepare_compartment_column(adata, source_col=location_col)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.target:
        passes = [{"title_format": args.title_format, "title_size": args.title_size,
                   "col_gap": args.col_gap, "targets": args.target}]
    else:
        passes = [{**p, "targets": [parse_target(t) for t in p["targets"]]}
                  for p in MANUSCRIPT_PASSES]

    n_ok = 0
    failures: List[str] = []
    for run in passes:
        for target in run["targets"]:
            requested_ct = str(target["cancer_type"])
            sample_id = str(target["sample_id"])
            pemt_vmax = float(target["pemt_vmax"]) if target["pemt_vmax"] is not None else float(args.pemt_vmax)
            spot_size = float(target["spot_size"]) if target["spot_size"] is not None else float(args.spot_size)
            value_gamma = float(target["value_gamma"]) if target["value_gamma"] is not None else float(args.value_gamma)
            try:
                cancer_type = resolve_cancer_type(adata, args.cancer_col, requested_ct)
                check_sample_in_cancer_type(adata, args.cancer_col, cancer_type, slide_col, sample_id)
                out_pdf = outdir / safe_name(cancer_type) / f"{safe_name(sample_id)}__pemt_overlap_panels.pdf"
                slide = select_slide(adata, slide_value=sample_id, slide_col=slide_col)
                plot_sample_panel(
                    slide=slide,
                    sample_id=sample_id,
                    cancer_type=cancer_type,
                    out_pdf=out_pdf,
                    overlap_specs=overlap_specs,
                    pemt_col=args.pemt_col,
                    pemt_vmax=pemt_vmax,
                    img_key=args.img_key,
                    spot_size=spot_size,
                    alpha_img=args.alpha_img,
                    alpha_spot=args.alpha_spot,
                    cmap=args.cmap,
                    dpi=args.dpi,
                    panel_title_size=args.panel_title_size,
                    colorbar_tick_size=args.colorbar_tick_size,
                    legend_fontsize=args.legend_fontsize,
                    transform=args.transform,
                    clip_percentile=args.clip_percentile,
                    layer=args.layer,
                    use_raw=args.use_raw,
                    overlap_gamma=args.overlap_gamma,
                    value_gamma=value_gamma,
                    pemt_label=args.pemt_label,
                    pemt_label_size=args.pemt_label_size,
                    crop_to_spots=args.crop_to_spots,
                    crop_pad_frac=args.crop_pad_frac,
                    title_gap=args.title_gap,
                    title_size=run["title_size"],
                    title_label=cancer_labels.get(cancer_type) or cancer_labels.get(requested_ct),
                    title_format=run["title_format"],
                    row_gap=args.row_gap,
                    panel_height=args.panel_height,
                    col_gap=run["col_gap"],
                )
                n_ok += 1
                print(f"[OK]   {cancer_type} | {sample_id} | pEMT vmax={pemt_vmax:g} | spot={spot_size:g} "
                      f"| value_gamma={value_gamma:g} -> {out_pdf}", flush=True)
            except Exception as e:  # noqa: BLE001
                failures.append(f"{requested_ct} | {sample_id}: {e}")
                print(f"[FAIL] {requested_ct} | {sample_id}: {e}", flush=True)

    print(f"Done. {n_ok} slide(s) saved, {len(failures)} failed.", flush=True)
    if failures:
        for line in failures:
            print("  " + line, flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
