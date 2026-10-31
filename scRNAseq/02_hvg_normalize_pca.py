import scanpy as sc
import numpy as np
import pandas as pd
sc.settings.n_jobs=32

# AnnData 객체 불러오기
print("Loading AnnData object...")
adata = sc.read_h5ad('/path/to/h5ad_raw/concat_data/01_QC_concat.h5ad')
adata.obs_names_make_unique()

print("Preserving raw counts in 인덱스 중복 확인...")
duplicates = adata.obs_names.duplicated()
print(f"중복된 인덱스가 있는지 확인: {duplicates.sum()}개")

############################################

print("Preserving raw counts in adata.layers['counts']...")
adata.layers["counts"] = adata.X.copy()

print("Normalizing and log-transforming the data...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print("Identifying highly variable genes...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='Sample')

print(f"Number of highly variable genes: {np.sum(adata.var['highly_variable'])}")

adata.raw = adata

print("Filtering the data to include only highly variable genes...")
adata = adata[:, adata.var.highly_variable]

print("Scaling the data...")
sc.pp.scale(adata, max_value=10)

print("Running PCA...")
sc.tl.pca(adata, svd_solver='arpack')

print("Computing neighbors and running UMAP...")
sc.pp.neighbors(adata)
sc.tl.umap(adata) # UMAP(non batch corr ver)

print("Saving the processed AnnData object...")
adata.write_h5ad('/path/to/h5ad_raw/concat_data/02_hvg_pca_concat.h5ad')

print("All processes completed successfully.")
