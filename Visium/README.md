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
