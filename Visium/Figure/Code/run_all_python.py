#!/usr/bin/env python3
"""Runs the Python panels (all of Figure 3, Fig. 4(a,b,d,f,g), Fig. 5d, Ext. Fig. 3/5b/6/7/8/10 and Supp. Fig. 11-13)
in order, using only the officially distributed Source Data in Paper_info/Source_Data as input. Fig. 4c, Fig. 5(a,b) and
Ext. Fig. 1/2/4 are reproduced in R and handled by run_all.R / run_all_extended.R (for overlapping panels this script
takes precedence, i.e. the Python code taken from the Fig345 Source Data upload codebase (final Ver2)).

Individual scripts can also be run on their own, e.g. `python Code/Figure3/*.py`, and do not depend on each other
(exception: Fig4e_ExtDataFig9_slide_panels.py / ExtDataFig9_grid_3x3.py; the latter uses the output of the former).

Fig. 4e/Ext. Fig. 9 and Fig. 5c are not drawn from the officially distributed Source Data but need the original Visium h5ad.
They are run only when the PCASSO_VISIUM_H5AD environment variable is set, and skipped otherwise.
"""
import os
import subprocess
import sys
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent

SOURCE_DATA_SCRIPTS = [
    "Figure3/Fig3d_ExtDataFig6_8_networks.py",
    "Figure3/Fig3ef_ExtDataFig7ef_centrality.py",
    "Figure3/Fig3cgh_ExtDataFig7cgh.py",
    "Figure3/Fig3b_ExtDataFig7b.py",
    "Figure3/ExtDataFig5b_7d_pair_categories.py",
    "Figure3/Suppfig11_bc_by_compartment.py",
    "Figure4/Fig4ab_signature_panels.py",
    "Figure4/Fig4d_ranked_scatter.py",
    "Figure4/Fig4fg_SupFig12-13_survival.py",
    "Figure4/ExtDataFig10_bulk_correlation.py",
    "Figure5/Fig5d_sc_vs_visium_scatter.py",
    "ExtendedFigure3.py",
]

H5AD_SCRIPTS = [
    "Figure4/Fig4e_ExtDataFig9_slide_panels.py",
    "Figure4/ExtDataFig9_grid_3x3.py",
    "Figure5/Fig5c_spatial_lr_overlap.py",
]


def run(rel_path: str) -> None:
    print(f"==== Running {rel_path} ====")
    subprocess.run([sys.executable, str(CODE_DIR / rel_path)], check=True, cwd=CODE_DIR)


if __name__ == "__main__":
    for rel_path in SOURCE_DATA_SCRIPTS:
        run(rel_path)

    if os.environ.get("PCASSO_VISIUM_H5AD"):
        for rel_path in H5AD_SCRIPTS:
            run(rel_path)
    else:
        print("PCASSO_VISIUM_H5AD not set: skipping Fig.4e/ExtFig.9/Fig.5c "
              f"({', '.join(H5AD_SCRIPTS)}) — these need the deposited Visium h5ad, not Source Data.")
