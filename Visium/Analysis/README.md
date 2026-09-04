# Analysis

Scripts that generate the intermediate per-slide and per-cell-type results consumed
by `Visium/Figure/`.

```
Analysis/
  Run_Cell2Location/
    make_singlecell_reference.py  — trains the scRNA-seq reference signature (NB
                                     regression) used by cell2location.
    run_cell2location.py          — per-organ Visium deconvolution against that
                                     reference; computes major-lineage (*_enriched)
                                     proportions.
  Run_STopover/
    run_stopover.py               — STopover ligand-receptor topological-similarity
                                     (Jaccard) analysis per slide, whole-slide or
                                     per compartment (Malignant/Boundary/Normal).
    calc_ct_lr_jaccard.R          — combines STopover LR results with cell2location
                                     cell-type connected components to compute
                                     cell-type-pair x LR-pair Jaccard indices.
```

All paths in these scripts are example directory names (`/path/to/...`) — edit them
for your environment before running. Subtype names use the post-curation cell type
labels (see `Data_Info/celltype_mapping.csv` in the analysis project for the
before/after mapping); Esophagus (ESCA) and Bladder are excluded, matching the
manuscript's 13-cancer-type (nonESCA) cohort.
