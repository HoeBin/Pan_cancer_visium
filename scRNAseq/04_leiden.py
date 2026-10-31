import platform
import numpy as np
import pandas as pd
import scanpy as sc

print("Loading AnnData object...")
adata_z=sc.read_h5ad('/path/to/h5ad_raw/concat_data/04_batch_corr_Harmony_clus.h5ad')
print("AnnData object loaded successfully.")
print(f"세포 개수: {len(adata_z.obs)}개")

doublet=pd.read_csv('/path/to/qc_output/scrublet_results_new.csv',index_col=0)
adata_z.obs['scrublet_score'] = doublet['scrublet_doublet_scores'].copy()
adata_z.obs['scrublet_predict'] = doublet['scrublet_predicted_doublets'].copy()
print(f"load doublet result.")

adata_z=adata_z[adata_z.obs['scrublet_score'] < 0.2]
print(f"더블렛 제거 후 세포 개수: {len(adata_z.obs)}개")

adata_z = adata_z[:, adata_z.var.highly_variable]

res = 0.2
print(res)
igraph = 'igraph'
print(igraph)
n_iterations=2
print(n_iterations)

sc.tl.leiden(adata_z, resolution=res, key_added=f'leiden_{res}',flavor=igraph, n_iterations=n_iterations)
print("Leiden clustering completed.")

print("Saving the AnnData object...")
adata_z.write_h5ad('/path/to/h5ad_raw/concat_data/06_batch_corr_Harmony_clus_lei.h5ad')
print("AnnData object saved successfully.")