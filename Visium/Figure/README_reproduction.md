# Figure reproduction package (Main Figures 1-5 + Extended Data Figures 1-10 + Supplementary Fig. 11-13)

This repository holds the code (`Code/`) and its output (`Output/`) that reproduce the data-driven
panels of the Main/Extended Data Figures using only the officially distributed
`Paper_info/Paper_info/Source_Data/*.xlsx` as input. R and Python are mixed — each panel keeps the
language of the codebase in which it was originally verified (see "Mixing R and Python" below).
`Paper_info/Paper_info` stores the original manuscript distribution material (Source Data, Manuscript,
Figure PDFs, etc.) as is, and this code writes no files into it.

Of the two code layers, only the **Figure layer** is included: there are only scripts that read the Source
Data (xlsx) and draw the figures, and the Analysis layer that recomputes statistics from external raw data
such as raw h5ad/TCGA counts is not included — this repository is a self-contained reproduction package
that uses only the official Source Data as input. The only exceptions are the three image panels that can be
drawn only with an h5ad file (Fig.4e, Ext.Fig.9, Fig.5c); these are included, and the need for an external
Visium h5ad path at run time is noted below.

## Directory structure

```
Visium/Figure/
  Code/
    common/config.R, config.py        — paths, palettes, Source Data readers (R/Python pair)
    Data_Info/subtype_lineage_map.csv — for the network node colors of Fig.3d/Ext.Fig.6/8 (93 subtypes -> 7 lineages)
    Figure1.R, Figure2.R                       — R, standalone
    Figure3/Fig3a_ExtDataFig7a_correlation_heatmap.R — R, Fig.3a + Ext.Fig.7a (Fig.3a uses Input/fig3a_source_data.csv, see "Fig.3a input" below)
    Figure3/Input/fig3a_source_data.csv        — bundled input: Fig.3a correlation matrix (not Source Data)
    Figure3/Fig3d_ExtDataFig6_8_networks.py, Fig3ef_ExtDataFig7ef_centrality.py,
           Fig3cgh_ExtDataFig7cgh.py, Fig3b_ExtDataFig7b.py,
           ExtDataFig5b_7d_pair_categories.py, Suppfig11_bc_by_compartment.py
                                                — Python, Fig.3(b-h) + Ext.Fig.5b/6/7/8 + Supp.Fig.11
    Figure4/Fig4ab_signature_panels.py, Fig4d_ranked_scatter.py,
           Fig4fg_SupFig12-13_survival.py, ExtDataFig10_bulk_correlation.py
                                                — Python, Fig.4(a,b,d), Fig.4(f,g)+Supp.Fig.12-13, Ext.Fig.10
    Figure4/Fig4c_pEMT_abundance_correlation.R — R, Fig.4c
    Figure4/Fig4e_ExtDataFig9_slide_panels.py,
           ExtDataFig9_grid_3x3.py  ⚠ needs an external h5ad; the latter uses the output of the former
                                                — Fig.4e, Ext.Fig.9
    Figure4/Input/pemt_scores.csv.gz            — bundled input used by the h5ad scripts above (not Source Data)
    Figure5/Fig5ab_pancancer_LR_zscore.R        — R, Fig.5a-b
    Figure5/Fig5d_sc_vs_visium_scatter.py       — Python, Fig.5d
    Figure5/Fig5c_spatial_lr_overlap.py  ⚠ needs an external h5ad — Fig.5c
    Figure5/Input/selected_slides.csv           — bundled input used by the h5ad script above (not Source Data)
    ExtendedFigure1.R, ExtendedFigure2.R, ExtendedFigure4.R, ExtendedFigure5.R (5a only) — R
    ExtendedFigure3.py                          — Python, companion of Fig.2 (Ext.Fig.3a-h)
    run_all.R            — runs the R panels in order (Fig1, 2, 3a, 4c, 5a-b)
    run_all_extended.R   — runs the R Extended panels in order (ExtFig 1, 2, 4, 5a)
    run_all_python.py    — runs the Python panels in order (all of Fig.3, Fig.4(a,b,d,f,g), Fig.5d, Ext.Fig.3/5b/6/7/8/10,
                            Supp.11-13; also continues with the h5ad panels when PCASSO_VISIUM_H5AD is set)
  Output/FigureN/{Plots,Tables}/, ExtendedFigureN/{Plots,Tables}/, SupplementaryFigureN/{Plots,Tables}/
  Paper_info/Paper_info/         — officially distributed manuscript material (input, read-only)
    Source_Data/*.xlsx           — only Source_Data is included in this repository (the manuscript, Figure PDFs/PPTX, etc. are excluded)
```

## Mixing R and Python, and precedence of overlapping panels

Fig.3/4/5 were built from two original codebases that partly overlap (the R-based `Github_code` and the
Python-based Fig345 Source Data upload codebase, final Ver2). For the overlapping panels, the Python code that
was verified against the original pipeline down to the pixel level on 2026-09-07 takes precedence:

| Panel | Adopted | Note |
|---|---|---|
| Fig.4a, 4b | Python (`Figure4/Fig4ab_signature_panels.py`) | The same 268 slides are matched across the 3 compartments (paired data), so a paired Wilcoxon test (+ Friedman) is used. The earlier R version had no `paired=TRUE` in `stat_compare_means` and therefore computed an unpaired Wilcoxon test. |
| Ext.Fig.5b | Python (`Figure3/ExtDataFig5b_7d_pair_categories.py`) | Paired Wilcoxon for the same reason. The R version (Ext.Fig.5b) was unpaired in the same way. |
| Ext.Fig.6 (Malignant/Normal co-enrichment network) | Python (`Figure3/Fig3d_ExtDataFig6_8_networks.py`) | The same script as Fig.3d/Ext.Fig.8 reproduces the node layout identically to the original, matched to Excel's 15-digit precision. |

Panels that do not overlap and keep their original code: Fig.4c (R, `Figure4/Fig4c_pEMT_abundance_correlation.R`), Fig.5a-b (R,
`Figure5/Fig5ab_pancancer_LR_zscore.R`), Ext.Fig.5a (R, `ExtendedFigure5.R`), Fig.1/2 and Ext.Fig.1/2/4 (R, none of which
are covered by the Python codebase).

## Scope and limitations (panels excluded because no Source Data exists)

- Figure 1b (BioRender workflow schematic)
- Figure 2a-d, 2f (scRNA-seq/Visium UMAP embedding, marker dot plot, malignancy UMAP density)
- Fig2g is reproduced only as an approximate dot plot because of the limited information in the Source Data (top 10 genes per compartment, dot size = |logFC|, color = score).
- Extended Fig.1a-b, 1d-e (scRNA-seq UMAP, marker gene dot plot); only 1c/1f (Hallmark pathway heatmap) are reproduced.
- Extended Fig.6 is not the manual Cytoscape layout but a force-directed layout recomputed from the edge list in the Source Data, so node coordinates may differ, but the layout is reproduced identically to the original at Excel's 15-digit precision (the network structure is identical).

## Running

Python environment: numpy, pandas, scipy, matplotlib, networkx, seaborn, openpyxl, h5py, lifelines (Fig.4f/g),
adjustText (Ext.Fig.10, Ext.Fig.9), pypdf (Ext.Fig.9 3x3 grid), anndata/scanpy (h5ad panels only). Arial is used
if it is installed, and the default sans-serif is used otherwise.

```bash
cd Visium/Figure
Rscript Code/run_all.R                    # Fig.1, 2, 3a (+Ext.Fig.7a), 4c, 5a-b
Rscript Code/run_all_extended.R           # Ext.Fig.1, 2, 4, 5a
python Code/run_all_python.py             # all of Fig.3, Fig.4(a,b,d,f,g), Fig.5d, Ext.Fig.3/5b/6/7/8/10, Supp.11-13

# Three image panels that can be drawn only with an h5ad (Fig.4e, Ext.Fig.9, Fig.5c) — optional
export PCASSO_VISIUM_H5AD=/real/path/_visium_ALL_260730.h5ad
python Code/run_all_python.py             # with the environment variable above set, this also runs the h5ad panels
```

## Fig.3a input

The `Fig.3a` sheet of `Source_Data_Fig3.xlsx` holds the old proportion-based correlation values (-0.62 to 0.61), which do not
match the manuscript text ("Epithelial cells showed near-zero correlations"). Fig.3a therefore uses
`Code/Figure3/Input/fig3a_source_data.csv`, which holds the abundance-based correlation (13 cancer types without ESCA,
0.016 to 0.833), as its default input. Once the Excel file is updated,

```bash
Rscript Code/Figure3/Fig3a_ExtDataFig7a_correlation_heatmap.R --fig3a-source=xlsx
```

uses the `Fig.3a` sheet directly. In the default (csv) mode, the maximum absolute difference between the CSV and the xlsx
sheet is printed on every run, so you can tell whether the sheet has been updated. Ext.Fig.7a reads the `ED Fig.7a` sheet of
`Source_Data_Extended_Fig7.xlsx` (BRCA downsampled to 30 slides, already abundance-based) as is.

Each script can also be run on its own (`Rscript Code/Figure1.R`, `python Code/Figure3/Fig3d_ExtDataFig6_8_networks.py`, etc.).
`Figure4/ExtDataFig9_grid_3x3.py` uses the output of `Figure4/Fig4e_ExtDataFig9_slide_panels.py`, so that
script has to be run first.
