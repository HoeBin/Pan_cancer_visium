# Visium

Code for processing, deconvolving, and analyzing the pan-cancer Visium spatial transcriptomics cohort in P-CASSO.

Corresponds to the manuscript Methods sections:
- **Visium data processing** — spot-level QC (≥500 UMIs, <30% mito), slide-level QC,
  joint embedding across FF/FFPE sections (BBKNN, UMAP).
- **Visium data deconvolution** — reference-based deconvolution with cell2location against
  the shared scRNA-seq reference, run per cancer type.
- **CNV-based tumor region annotation** — Malignant/Boundary/Normal compartment
  segmentation with Cottrazm (SME normalization + inferCNV).
- **Spatial co-localization** — Jaccard-based co-localization index (STopover) at the
  whole-slide and compartment level.
- **Spatial autocorrelation analysis** — Moran's I per cell type per slide.
- **Cell-type co-enrichment network analysis** — compartment-specific networks and
  betweenness centrality (Cytoscape).
- **Spatial ligand-receptor association** — co-localization of curated LR pairs
  (CellTalkDB) with cell-type pairs.
- **Gene-set scoring of Visium spots** — epithelial, mesenchymal, and pEMT signature
  scoring.

## Figure reproduction

[`Figure/`](Figure/README_reproduction.md) reproduces the data-driven panels of the main and extended figures
(R and Python) from the officially distributed Source Data (`Figure/Paper_info/Paper_info/Source_Data/*.xlsx`).
Outputs (CSV tables and PDF plots) are written to `Figure/Output/`.

```bash
cd Visium/Figure
Rscript Code/run_all.R                    # Fig. 1, 2, 3a (+ Ext. Fig. 7a), 4c, 5a-b
Rscript Code/run_all_extended.R           # Ext. Fig. 1, 2, 4, 5a
python Code/run_all_python.py             # Fig. 3, 4, 5d and the remaining Ext./Supp. Figures
```

Fig. 3a uses the bundled `Figure/Code/Figure3/Input/fig3a_source_data.csv` (abundance-based correlation,
13 cancer types); see [`Figure/README_reproduction.md`](Figure/README_reproduction.md) for the input details,
the panels that are not reproducible from Source Data, and the optional h5ad-based panels.
