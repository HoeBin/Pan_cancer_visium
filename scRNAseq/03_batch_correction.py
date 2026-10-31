import scanpy as sc
import pandas as pd
import numpy as np
import bbknn
sc.settings.n_jobs=32
import os, sys


print("Loading AnnData object...")
adata = sc.read_h5ad('/path/to/h5ad_raw/concat_data/02_hvg_pca_concat.h5ad')
duplicates = adata.obs_names.duplicated()
print(f"중복된 인덱스가 있는지 확인: {duplicates.sum()}개")

print("load log_normalizeed expression...")
adata=adata.raw.to_adata()
adata.raw = adata

print("Running PCA...")
sc.tl.pca(adata)

print("Batch correction...")
sc.external.pp.harmony_integrate(adata, key='batch_key', max_iter_harmony=50) # batch_key : Chemistry
bbknn.bbknn(adata, batch_key='Sample',use_rep='X_pca_harmony',n_pcs=40)

print("Calculating UMAP...")
sc.tl.umap(adata)  # UMAP(batch corrrection ver)
print("UMAP calculation completed.")

# save
print("Saving the results...")

adata.write_h5ad('/path/to/h5ad_raw/concat_data/04_batch_corr_Harmony_clus.h5ad')
print("All processes completed successfully.")