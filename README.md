# P-CASSO: A compartment-resolved pan-cancer spatial atlas of intercellular interactions and prognostic ligand-receptor axes

This repository contains the analysis code accompanying the manuscript:

> **P-CASSO: A compartment-resolved pan-cancer spatial atlas of intercellular interactions and prognostic ligand-receptor axes**
> Jaewoo Mo, Hoebin Chung, Gahyun Kim, Gyeong-Jin Shin, Jeongbin Park, Woong-Yang Park, Junil Kim
> Web portal: [https://p-casso.ai](https://p-casso.ai)

## Overview

P-CASSO is a compartment-resolved pan-cancer spatial atlas that pairs a single-cell RNA-seq
reference (686 samples; 2,183,305 cells; 93 fine-grained cell types) with 10x Visium data
(268 slides; 635,720 spots) across 13 cancer types. Every Visium slide is deconvolved against
the shared reference and segmented into **Malignant**, **Boundary**, and **Normal**
compartments using morphology- and CNV-informed scoring. Cell-type co-enrichment is
quantified with a Jaccard-based spatial co-localization index, used to build compartment-level
co-enrichment networks and to prioritize ligand-receptor (LR) interactions, whose prognostic
value is further tested in independent TCGA cohorts.

## Analysis pipeline

The code in this repository implements the following stages, corresponding to the Methods
section of the manuscript:

1. **scRNA-seq reference construction** — QC, doublet removal (Scrublet), integration and
   batch correction (Harmony), lineage/subtype clustering (Leiden) and marker-based
   annotation into 93 fine-grained cell types.
2. **Visium processing** — spot/slide-level QC, joint embedding (BBKNN, UMAP) across
   fresh-frozen and FFPE sections.
3. **Visium deconvolution** — reference-based deconvolution with cell2location, run
   per cancer type against the shared scRNA-seq reference.
4. **CNV-based compartment annotation** — Malignant/Boundary/Normal segmentation with
   Cottrazm (SME normalization + inferCNV).
5. **Spatial co-localization & network analysis** — Jaccard-based co-localization index
   (STopover), compartment-specific co-enrichment networks and betweenness centrality
   (Cytoscape).
6. **Ligand-receptor prioritization** — spatial co-localization of curated LR pairs
   (CellTalkDB) within compartment- and cell-type-pair-defined neighborhoods.
7. **Gene-set scoring** — epithelial, mesenchymal and pEMT signature scoring of Visium
   spots (Scanpy `score_genes`).
8. **TCGA bulk deconvolution & survival analysis** — BayesPrism/OmicVerse deconvolution
   of TCGA bulk RNA-seq, Cox proportional-hazards modeling and Kaplan-Meier analysis
   (lifelines) of LR/cell-abundance signals across cancer types.

## Key dependencies

- Python: `scanpy`, `cell2location`, `harmonypy`, `bbknn`, `scrublet`, `Cottrazm`,
  `infercnvpy`/`inferCNV`, `STopover`, `gseapy`, `omicverse`, `lifelines`, `statsmodels`
- R: `TCGAbiolinks`, `spdep`
- Cytoscape (network visualization and betweenness centrality)

## Data availability

- scRNA-seq and Visium datasets analyzed in this study are publicly available (GEO, 10x
  Genomics, Zenodo, EMBL-EBI, NGDC/GSA, Mendeley Data, Single Cell Portal); accession
  numbers are listed in the manuscript's Supplementary Tables 1–4.
- TCGA bulk RNA-seq and clinical data were obtained from the NCI Genomic Data Commons
  (GDC).
- The integrated pan-cancer single-cell reference and processed spatial outputs are available
  through the [P-CASSO web portal](https://p-casso.ai) and Zenodo.

## Citation

If you use this code or the P-CASSO atlas, please cite the manuscript above.
