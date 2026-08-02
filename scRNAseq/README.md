# scRNA-seq

Code for constructing the pan-cancer single-cell RNA-seq reference used throughout P-CASSO.

Corresponds to the manuscript Methods sections:
- **Single-cell RNA-seq reference construction** — QC (≥1,000 UMIs, ≥500 genes, <20% mito,
  ≤7,000 genes), doublet removal (Scrublet), normalization, HVG selection, PCA, batch
  correction (Harmony), Leiden clustering, and major-lineage annotation (epithelial, T cell,
  B cell, myeloid, endothelial, fibroblast, mural).
- Per-lineage subclustering (Harmony + BBKNN, Leiden) and fine-grained annotation into 93
  cell types.
- **Cell-type composition estimates** — per-sample lineage fractions used for the
  scRNA-seq vs. Visium composition comparison.
- **Pathway over-representation analysis** — marker gene over-representation against
  MSigDB Hallmark gene sets (GSEApy/Enrichr).
