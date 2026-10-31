import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scrublet as scr
import session_info
import bbknn

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=200, facecolor='white', color_map = 'RdPu',)

#loading adata
adata_endo=sc.read_h5ad(
    '/path/to/h5ad_raw/concat_data/07_fibro+mural.h5ad')

sample_counts = adata_endo.obs['Sample'].value_counts()
samples_to_remove = sample_counts[sample_counts < 10].index
adata_endo = adata_endo[~adata_endo.obs['Sample'].isin(samples_to_remove)].copy()
sc.pp.highly_variable_genes(adata_endo, 
                            min_mean=0.0125, max_mean=3, min_disp=0.5,
                            n_top_genes=3000,
                            batch_key='Sample'
                           )
print(sum(adata_endo.var.highly_variable))
sc.pl.highly_variable_genes(adata_endo)
print(f"detected {np.sum(adata_endo.var['highly_variable'])} highly variable genes")
adata_endo.raw = adata_endo
adata_endo = adata_endo[:, adata_endo.var.highly_variable]
sc.pp.scale(adata_endo, max_value=10)
sc.tl.pca(adata_endo)
sc.external.pp.harmony_integrate(adata_endo, key = 'Sample', max_iter_harmony = 50)
bbknn.bbknn(adata_endo, batch_key='batch_key',use_rep='X_pca_harmony',n_pcs=40)
sc.tl.umap(adata_endo)
adata_endo.write_h5ad('/path/to/h5ad_subcluster/07_fibro+mural_sample_batchkey_hvg3000.h5ad')