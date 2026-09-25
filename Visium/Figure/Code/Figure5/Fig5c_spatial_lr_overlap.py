#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fig. 5c: spatial co-enrichment (two cell types) and co-expression (ligand, receptor) of cancer-type-specific
LR pairs on five Visium slides, drawn over the hires H&E image.

Input : deposited Visium AnnData _visium_ALL_260730.h5ad (env PCASSO_VISIUM_H5AD or --h5ad; not bundled) and
        the bundled manifest Input/Fig5c/selected_slides.csv (--manifest: sample_id, cancer_type, cell_type_1,
        cell_type_2, ligand, receptor, jaccard). Cell-type abundance = obs "<cell type>_enriched" x obs "sum";
        expression from X; coordinates = obsm "spatial" x tissue_hires_scalef of uns["spatial"][sample_id].
Output: Output/Figure5/Plots/fig5c_<cancer_type>_<sample_id>_<ligand>-<receptor>.pdf,
        one per slide (five in total).
Each feature is arcsinh-transformed, clipped at the 99th percentile and scaled to 0-1; spots mix red
(first feature only), blue (second only) and green (overlap, the smaller value), brightness = the larger
value, with a 2-D colour key per panel. The jaccard column is shown in the title only. Text uses the
matplotlib default font (no Arial setup).
"""
from __future__ import annotations

# "Libraries and paths" ----
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # Code/ for `common`
from common import config as cfg

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset


# "Inputs" ----
INPUT_DIR = cfg.FIG5C_INPUT                      # bundled: selected_slides.csv
MANIFEST_CSV = INPUT_DIR / "selected_slides.csv"
VISIUM_H5AD = cfg.VISIUM_H5AD                    # external: PCASSO_VISIUM_H5AD env var or --h5ad

# "Outputs" ----
# fig5c_<cancer_type>_<sample_id>_<ligand>-<receptor>.pdf -> cfg.output_dirs("Fig5")[0]
# (Output/Figure5/Plots/; resolved inside main() so importing/--help creates nothing).

CLIP_PERCENTILE = 99.0
SPOT_SIZE = 9.0
FIGSIZE = (13.8, 5.4)
DPI = 700
KEY_GRID = 80

COLOR_LEFT = np.array([0.90, 0.12, 0.14])
COLOR_OVERLAP = np.array([0.10, 0.75, 0.20])
COLOR_RIGHT = np.array([0.10, 0.35, 0.95])


def scale_feature(v: np.ndarray) -> np.ndarray:
    out = np.arcsinh(np.nan_to_num(np.asarray(v, dtype=float).ravel(), nan=0.0, posinf=0.0, neginf=0.0))
    hi = np.nanpercentile(out, CLIP_PERCENTILE)
    if hi <= 0:
        return np.zeros_like(out)
    return np.clip(out, 0, hi) / hi


def mix_rgb(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left = np.clip(np.asarray(left, dtype=float).ravel(), 0.0, 1.0)
    right = np.clip(np.asarray(right, dtype=float).ravel(), 0.0, 1.0)
    w_green = np.minimum(left, right)
    w_red = np.clip(left - w_green, 0.0, 1.0)
    w_blue = np.clip(right - w_green, 0.0, 1.0)
    w_sum = w_red + w_green + w_blue

    rgb = np.zeros((left.size, 3), dtype=float)
    valid = w_sum > 1e-12
    rgb[valid] = (
        w_red[valid, None] * COLOR_LEFT
        + w_green[valid, None] * COLOR_OVERLAP
        + w_blue[valid, None] * COLOR_RIGHT
    ) / w_sum[valid, None]
    rgb = rgb * np.maximum(left, right)[:, None]
    return np.clip(rgb, 0.0, 1.0)


def color_key() -> np.ndarray:
    grid = np.linspace(0.0, 1.0, KEY_GRID)
    left, right = np.meshgrid(grid, grid)
    return mix_rgb(left.ravel(), right.ravel()).reshape(KEY_GRID, KEY_GRID, 3)


def load_slide(h5ad: Path, sample_id: str, cell_types, genes) -> tuple[np.ndarray, np.ndarray, dict]:
    """Return hires-pixel spot coordinates, the hires image and per-spot feature vectors."""
    with h5py.File(h5ad, "r") as f:
        rows = np.flatnonzero(np.asarray(read_elem(f["obs"]["sample_id"])) == sample_id)
        total = read_elem(f["obs"]["sum"])[rows]
        values = {ct: read_elem(f["obs"][f"{ct}_enriched"])[rows] * total for ct in cell_types}
        var_index = read_elem(f["var"]).index
        x = sparse_dataset(f["X"])[rows]
        values.update({g: x[:, var_index.get_loc(g)].toarray().ravel() for g in genes})
        spatial = read_elem(f["obsm"]["spatial"])[rows]
        lib = read_elem(f["uns"]["spatial"][sample_id])
    xy = spatial * lib["scalefactors"]["tissue_hires_scalef"]
    return xy, lib["images"]["hires"], values


def draw_slide(xy: np.ndarray, img: np.ndarray, values: dict, pairs, title: str, out_pdf: Path) -> None:
    """pairs: ((cell_type_1, cell_type_2), (ligand, receptor)); names are keys of values."""
    key = color_key()

    fig = plt.figure(figsize=FIGSIZE)
    gs = fig.add_gridspec(1, 7, width_ratios=[1.0, 0.10, 0.23, 0.30, 1.0, 0.10, 0.23], wspace=0.0)
    slide_axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 4])]
    key_axes = [fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[0, 6])]
    panels = zip(slide_axes, key_axes, pairs, ["group_pair", "lr_pair"], ["Co-enrichment", "Co-expression"])

    for ax, cax, (left, right), prefix, key_title in panels:
        rgb = mix_rgb(scale_feature(values[left]), scale_feature(values[right]))
        ax.imshow(img)
        ax.scatter(xy[:, 0], xy[:, 1], s=SPOT_SIZE, c=rgb, edgecolors="none")
        ax.set_axis_off()
        ax.set_title(f"{prefix}: {left} vs {right}", fontsize=11)

        cax.imshow(key, origin="lower", extent=[0, 1, 0, 1], interpolation="nearest")
        cax.set_aspect("equal")
        cax.set_title(key_title, fontsize=17, pad=14)
        cax.set_xlabel(left, fontsize=16, labelpad=10)
        cax.set_ylabel(right, fontsize=16, labelpad=10)
        cax.set_xticks([0.0, 0.5, 1.0])
        cax.set_yticks([0.0, 0.5, 1.0])
        cax.yaxis.set_label_position("right")
        cax.yaxis.tick_right()
        cax.tick_params(axis="both", labelsize=15, length=2, pad=4)
        for spine in cax.spines.values():
            spine.set_linewidth(0.8)

    fig.suptitle(title, fontsize=12)
    fig.subplots_adjust(left=0.03, right=0.985, bottom=0.08, top=0.90)
    fig.savefig(out_pdf, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Fig. 5c spatial co-enrichment / co-expression panels.")
    ap.add_argument("--h5ad", default=str(VISIUM_H5AD),
                    help="Visium AnnData _visium_ALL_260730.h5ad (not bundled; cfg.VISIUM_H5AD; default: %(default)s)")
    ap.add_argument("--manifest", default=str(MANIFEST_CSV),
                    help="selected_slides.csv manifest (cfg.FIG5C_INPUT; default: %(default)s)")
    args = ap.parse_args()

    h5ad = cfg.require_external(Path(args.h5ad), "Visium AnnData (_visium_ALL_260730.h5ad)")
    fig_dir = cfg.output_dirs("Fig5")[0]           # Output/Figure5/Plots/ (created here)
    manifest = pd.read_csv(args.manifest)
    for row in manifest.itertuples(index=False):
        pairs = ((row.cell_type_1, row.cell_type_2), (row.ligand, row.receptor))
        xy, img, values = load_slide(h5ad, row.sample_id, pairs[0], pairs[1])
        title = (
            f"cancer_type={row.cancer_type} | sample={row.sample_id} | "
            f"group_pair={row.cell_type_1}-{row.cell_type_2} | lr_pair={row.ligand}-{row.receptor} | "
            f"jaccard={row.jaccard:.4f}"
        )
        out_pdf = fig_dir / f"fig5c_{row.cancer_type}_{row.sample_id}_{row.ligand}-{row.receptor}.pdf"
        draw_slide(xy, img, values, pairs, title, out_pdf)
        print(out_pdf.name)


if __name__ == "__main__":
    main()
