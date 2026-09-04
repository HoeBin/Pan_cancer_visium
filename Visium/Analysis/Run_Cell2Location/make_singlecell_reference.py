"""Single-cell reference signature builder (NB regression) for cell2location.

Consolidates Make_SingleCell_Reference_subnum.ipynb ("Run_Cell2Location" version)
into one script. Produces inf_aver.csv, the reference signature matrix that
run_cell2location.py reads from REF_RUN_NAME. Run this once, before running
run_cell2location.py for any organ.

Usage:
    python make_singlecell_reference.py --gpu 2
"""

import argparse
import os

import cell2location
import pandas as pd
import scanpy as sc
import torch
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes
from matplotlib import rcParams

rcParams["pdf.fonttype"] = 42

# "Paths" -----------------------------------------------------------------
# The paths below are example directory names. Edit them for your environment
# before running.
PROJECT_DIR = "/path/to/project_root"
REF_ADATA_PATH = (
    "/path/to/input/scRNAseq_reference/"
    "annotated_reference.h5ad"
)
REF_RUN_NAME = (
    "/path/to/output/pan_visium_results/"
    "single_cell_reference/reference_signatures"
)
# consumed downstream (e.g. Check_SlideGens.ipynb) via this fixed absolute path
TOTAL_GENES_PATH = os.path.join(
    PROJECT_DIR, "project_data/scRNAseq_totalgenes.txt"
)

# "QC / model parameters" ---------------------------------------------------
GENE_FILTER_PARAMS = dict(
    cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12
)
TRAIN_PARAMS = dict(max_epochs=250, use_gpu=True)
EXPORT_QUANTILES = ["q05", "q50", "q95", "q0001"]


def load_reference():
    adata_ref = sc.read_h5ad(REF_ADATA_PATH)
    adata_ref.var["SYMBOL"] = adata_ref.var.index
    del adata_ref.raw
    return adata_ref


def filter_reference_genes(adata_ref):
    selected = filter_genes(adata_ref, **GENE_FILTER_PARAMS)
    return adata_ref[:, selected].copy()


def train_regression_model(adata_ref, gpu):
    cell2location.models.RegressionModel.setup_anndata(
        adata=adata_ref, batch_key="Sample", labels_key="sub_anno"
    )
    mod = RegressionModel(adata_ref)

    torch.cuda.set_device(gpu)
    mod.train(**TRAIN_PARAMS)

    adata_ref = mod.export_posterior(
        adata_ref,
        use_quantiles=True,
        add_to_varm=EXPORT_QUANTILES,
        sample_kwargs={"batch_size": 2500, "use_gpu": True},
    )
    return adata_ref, mod


def extract_inf_aver(adata_ref):
    factor_names = adata_ref.uns["mod"]["factor_names"]
    cols = [f"means_per_cluster_mu_fg_{i}" for i in factor_names]
    if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][cols].copy()
    else:
        inf_aver = adata_ref.var[cols].copy()
    inf_aver.columns = factor_names
    return inf_aver


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gpu", type=int, default=2, help="CUDA device index")
    args = parser.parse_args()

    # matches the fixed value used in the notebook (not tied to torch.cuda.set_device)
    os.environ["THEANO_FLAGS"] = "device=cuda0,floatX=float32,force_device=True"
    os.makedirs(REF_RUN_NAME, exist_ok=True)

    print("# Loading reference scRNA-seq data")
    adata_ref = load_reference()
    print(f"  {adata_ref.shape}")

    print("# Saving full reference gene list")
    pd.DataFrame(adata_ref.var_names.tolist()).to_csv(
        TOTAL_GENES_PATH, index=False, header=False
    )

    print("# QC-filtering reference genes")
    adata_ref = filter_reference_genes(adata_ref)
    print(f"  {adata_ref.shape}")

    print("# Training NB regression model")
    adata_ref, mod = train_regression_model(adata_ref, args.gpu)

    print("# Saving model and reference anndata")
    mod.save(REF_RUN_NAME, overwrite=True)
    adata_ref.write(os.path.join(REF_RUN_NAME, "sc.h5ad"))

    print("# Extracting and saving inf_aver (reference signature matrix)")
    inf_aver = extract_inf_aver(adata_ref)
    inf_aver.to_csv(os.path.join(REF_RUN_NAME, "inf_aver.csv"))
    print(f"  {inf_aver.shape} -> {REF_RUN_NAME}/inf_aver.csv")


if __name__ == "__main__":
    main()
