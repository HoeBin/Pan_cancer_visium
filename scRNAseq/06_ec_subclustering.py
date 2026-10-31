import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scrublet as scr
import bbknn
import session_info

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=200, facecolor='white', color_map = 'RdPu',)

#loading adata
adata_endo=sc.read_h5ad(
    '/path/to/h5ad_raw/concat_data/07_ec.h5ad')


sample_counts = adata_endo.obs['Sample'].value_counts()
samples_to_remove = sample_counts[sample_counts < 10].index
adata_endo = adata_endo[~adata_endo.obs['Sample'].isin(samples_to_remove)]
sc.pp.highly_variable_genes(adata_endo, 
                            min_mean=0.0125, max_mean=3, min_disp=0.5,
                            n_top_genes=4000,
                            batch_key='Sample'
                           )
print(sum(adata_endo.var.highly_variable))
sc.pl.highly_variable_genes(adata_endo)
print(f"detected {np.sum(adata_endo.var['highly_variable'])} highly variable genes")

adata_endo.raw = adata_endo
adata_endo = adata_endo[:, adata_endo.var.highly_variable]
sc.pp.scale(adata_endo, max_value=10)
sc.tl.pca(adata_endo)
bbknn.bbknn(adata_endo, batch_key='Sample')
sc.tl.umap(adata_endo)
sc.tl.leiden(adata_endo,resolution=1.0,flavor="igraph", n_iterations=2)
bbknn.ridge_regression(adata_endo, batch_key=['batch_key'], confounder_key=['leiden'])
sc.pp.pca(adata_endo)
bbknn.bbknn(adata_endo, batch_key='Sample')
sc.tl.umap(adata_endo)
adata_endo.write_h5ad('/path/to/h5ad_subcluster/07_ec_bbk.h5ad')