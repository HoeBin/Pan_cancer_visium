"""Cell2location spatial deconvolution — single-organ runner.

Consolidates the 14 per-organ notebooks in this folder (Cell2Location_<Organ>.ipynb,
"Run_Cell2Location" version) into one parameterized script. Only organs corresponding
to cancer types analyzed in the manuscript are included; Bladder is excluded (see
Data_Info.r: "Remove Bladder cancer"). Esophagus (ESCA) is also excluded — the
published figures use the nonESCA cohort (13 cancer types).

Model/training parameters match what was actually run per organ, transcribed from
each notebook's non-commented cells — except training batch_size and GPU device
index, which are unified to a single value across all organs (5000 and 0).

Usage:
    python run_cell2location.py --organ Breast
    python run_cell2location.py --list-organs
"""

import argparse
import os

import anndata
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from matplotlib import rcParams

rcParams["pdf.fonttype"] = 42

import cell2location
from cell2location.utils import select_slide

# "Paths" -----------------------------------------------------------------
# The paths below are example directory names. Edit them for your environment
# before running.
PROJECT_DIR = "/path/to/project_root"
INPUT_DIR = "/path/to/input/visium_raw_data/"
REF_RUN_NAME = (
    "/path/to/output/pan_visium_results/"
    "single_cell_reference/reference_signatures"
)
RUN_NAME = "/path/to/output/pan_visium_results/cell2location_deconvolution/"
SLIDE_INFO_PATH = os.path.join(
    PROJECT_DIR, "project_data/Visium_slide_info.txt"
)

# "Manuscript cancer types" -------------------------------------------------
# Organ -> cancer type abbreviation used in the paper (Data_Info.r / OrganToType.txt).
# Bladder and Esophagus(ESCA) are intentionally excluded — not part of the
# manuscript's analyzed (nonESCA) cancer types.
ORGAN_CANCER_TYPE = {
    "Breast": "BRCA",
    "Colon": "COCA",
    "Endometrium": "UECA",
    "HeadnNeck_Oral": "HNCA",
    "Kidney": "KICA",
    "Liver": "LICA",
    "Lung": "LUCA",
    "Ovary": "OVCA",
    "Pancreas": "PACA",
    "Prostate": "PRCA",
    "Skin": "SKCA",
    "Stomach": "STCA",
    "Thyroid": "THCA",
}

# "Run settings" -------------------------------------------------------------
# The original notebooks used a different mod.train(batch_size=...) and
# torch.cuda.set_device(...) per organ (job scheduling artifacts, not scientific
# parameters). Unified to a single value across all organs.
TRAIN_BATCH_SIZE = 5000
DEFAULT_GPU = 0

# "Visium dataset sources per organ" ---------------------------------------
# Each entry is either:
#   {"glob": "<relative_dir>/"}                       -> all files under INPUT_DIR/<dir>
#   {"files": [...]}                                    -> explicit relative file list
#   {"glob": "<relative_dir>/", "keep_identifiers": [..]} -> glob, then keep files whose
#                                                            name contains any identifier
ORGAN_DATASETS = {
    "Breast": [
        {"glob": "10x/Breast/raw_h5ad/"},
        {
            "files": [
                "Itai_GSE203612/raw_h5ad/GSM6177599_NYU_BRCA0_Vis_raw.h5ad",
                "Itai_GSE203612/raw_h5ad/GSM6177601_NYU_BRCA1_Vis_raw.h5ad",
                "Itai_GSE203612/raw_h5ad/GSM6177603_NYU_BRCA2_Vis_raw.h5ad",
            ]
        },
        {"glob": "Breast_GSE210616/raw_h5ad/"},
        {"glob": "Breast_4739739/raw_h5ad/"},
        {"glob": "Breast_GSE242311/raw_h5ad/"},
        {"glob": "Breast_GSE243022/raw_h5ad/"},
    ],
    "Colon": [
        {"glob": "CRC_GSE226997/raw_h5ad/"},
        {"glob": "10x/Large_Intestine/raw_h5ad/"},
        {
            "files": [
                "CRC_Liver_GSE225857/raw_h5ad/GSM7058756_raw.h5ad",
                "CRC_Liver_GSE225857/raw_h5ad/GSM7058757_raw.h5ad",
                "CRC_Liver_GSE225857/raw_h5ad/GSM7058758_raw.h5ad",
                "CRC_Liver_GSE225857/raw_h5ad/GSM7058759_raw.h5ad",
            ]
        },
        {
            "files": [
                "scCRLM_Atlas/raw_h5ad/ST-colon1_raw.h5ad",
                "scCRLM_Atlas/raw_h5ad/ST-colon2_raw.h5ad",
                "scCRLM_Atlas/raw_h5ad/ST-colon3_raw.h5ad",
                "scCRLM_Atlas/raw_h5ad/ST-colon4_raw.h5ad",
            ]
        },
        {"glob": "CRC_7760264/raw_h5ad/"},
        {"files": ["Colon_10x/raw_h5ad/Colon_10x_Visium_SampleP2_CRC_raw.h5ad"]},
    ],
    "Endometrium": [
        {"files": ["Itai_GSE203612/raw_h5ad/GSM6177623_NYU_UCEC3_Vis_raw.h5ad"]},
        {
            "files": [
                "Endometrial_GSE225690/raw_h5ad/01_034_C1d1_raw.h5ad",
                "Endometrial_GSE225690/raw_h5ad/01_034_C3d1_raw.h5ad",
                "Endometrial_GSE225690/raw_h5ad/01_039_C1d1_raw.h5ad",
                "Endometrial_GSE225690/raw_h5ad/01_039_C3d1_raw.h5ad",
            ]
        },
    ],
    "HeadnNeck_Oral": [
        {
            "files": [
                "HeadAndNeck_GSE181300/raw_h5ad/GSM5494475_HNSCC201125T04_raw.h5ad",
                "HeadAndNeck_GSE181300/raw_h5ad/GSM5494476_HNSCC201125T05_raw.h5ad",
                "HeadAndNeck_GSE181300/raw_h5ad/GSM5494477_HNSCC201125T07_raw.h5ad",
                "HeadAndNeck_GSE181300/raw_h5ad/GSM5494478_HNSCC201125T10_raw.h5ad",
                "HeadAndNeck_GSE181300/raw_h5ad/GSM5494479_P210325T1_raw.h5ad",
                "HeadAndNeck_GSE181300/raw_h5ad/GSM5494480_P210325T3_raw.h5ad",
                "HeadAndNeck_GSE181300/raw_h5ad/GSM5494481_P210325T5_raw.h5ad",
                "HeadAndNeck_GSE181300/raw_h5ad/GSM5494482_P210325T6_raw.h5ad",
            ]
        },
        {"glob": "Oral_GSE208253/raw_h5ad/"},
        {"glob": "Oral_GSE220978/raw_h5ad/"},
        {"glob": "Nasopharyngeal_GSE200310/raw_h5ad/"},
    ],
    "Kidney": [
        {"glob": "Kidney_GSE175540/raw_h5ad/"},
    ],
    "Liver": [
        {"files": ["Itai_GSE203612/raw_h5ad/GSM6177612_NYU_LIHC1_Vis_raw.h5ad"]},
        {"glob": "Liver_GSE238264/raw_h5ad/"},
        {
            "files": [
                "CRC_Liver_GSE225857/raw_h5ad/GSM7058760_raw.h5ad",
                "CRC_Liver_GSE225857/raw_h5ad/GSM7058761_raw.h5ad",
            ]
        },
        {"glob": "Liver_GSE217414/raw_h5ad/"},
        {
            "files": [
                "scCRLM_Atlas/raw_h5ad/ST-liver1_raw.h5ad",
                "scCRLM_Atlas/raw_h5ad/ST-liver2_raw.h5ad",
                "scCRLM_Atlas/raw_h5ad/ST-liver3_raw.h5ad",
                "scCRLM_Atlas/raw_h5ad/ST-liver4_raw.h5ad",
            ]
        },
        {
            "files": [
                "Liver_10119273/raw_h5ad/Tumor1_raw.h5ad",
                "Liver_10119273/raw_h5ad/Tumor2_raw.h5ad",
            ]
        },
        {
            "files": [
                "Liver_HRA000437/raw_h5ad/HCC-1T_raw.h5ad",
                "Liver_HRA000437/raw_h5ad/HCC-2T_raw.h5ad",
                "Liver_HRA000437/raw_h5ad/HCC-3T_raw.h5ad",
                "Liver_HRA000437/raw_h5ad/HCC-4T_raw.h5ad",
            ]
        },
        {
            "files": [
                "Liver_GSE245908/raw_h5ad/GSM7850822_CHC20_raw.h5ad",
                "Liver_GSE245908/raw_h5ad/GSM7850823_CHC23_raw.h5ad",
            ]
        },
    ],
    "Lung": [
        {
            "files": [
                "10x/Lung/raw_h5ad/CytAssist_11mm_FFPE_Human_Lung_Cancer_230214_raw.h5ad",
                "10x/Lung/raw_h5ad/CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_220714_raw.h5ad",
            ]
        },
        {
            "glob": "Lung_E_MTAB_13530/raw_h5ad/",
            "keep_identifiers": [
                "P10_T1", "P10_T2", "P10_T3", "P10_T4",
                "P11_T1", "P11_T2", "P11_T3", "P11_T4",
                "P15_T1", "P15_T2", "P16_T1", "P16_T2",
                "P17_T1", "P17_T2", "P19_T1", "P19_T2",
                "P24_T1", "P24_T2", "P25_T1", "P25_T2",
            ],
        },
    ],
    "Ovary": [
        {"glob": "10x/Ovary/raw_h5ad/"},
        {
            "files": [
                "Itai_GSE203612/raw_h5ad/GSM6177614_NYU_OVCA1_Vis_raw.h5ad",
                "Itai_GSE203612/raw_h5ad/GSM6177617_NYU_OVCA3_Vis_raw.h5ad",
            ]
        },
        {"glob": "Ovarian_GSE227019/raw_h5ad/"},
        {"glob": "Ovarian_GSE213699/raw_h5ad/"},
        {"glob": "Ovary_GSE211956/raw_h5ad/"},
        {"glob": "Ovary_GSE224335/raw_h5ad/"},
    ],
    "Pancreas": [
        {"files": ["Itai_GSE203612/raw_h5ad/GSM6177618_NYU_PDAC1_Vis_raw.h5ad"]},
        {
            "files": [
                "IPMN_GSE233293/raw_h5ad/GSM7421790_PDAC_1_raw.h5ad",
                "IPMN_GSE233293/raw_h5ad/GSM7421791_PDAC_2_raw.h5ad",
                "IPMN_GSE233293/raw_h5ad/GSM7421792_PDAC_3_raw.h5ad",
            ]
        },
        {"glob": "Pancreas_GSE211895/raw_h5ad/"},
        {"glob": "Pancreas_GSE235315/raw_h5ad/"},
    ],
    "Prostate": [
        {"glob": "10x/Prostate/raw_h5ad/"},
        {"glob": "Prostate_10.17632_4w6krnywhn.1/raw_h5ad/"},
        {"glob": "Prostate_GSE230282/raw_h5ad/"},
    ],
    "Skin": [
        {"glob": "Skin_GSE144239/raw_h5ad/"},
        {"glob": "Skin_10.17632_2bh5fchcv6.1/raw_h5ad/"},
        {
            "glob": "Skin_E_MTAB_13084/raw_h5ad/",
            "keep_identifiers": [
                "WSSKNKCLsp12140272", "WSSKNKCLsp12887269", "WSSKNKCLsp12887267",
                "WSSKNKCLsp12140270", "WSSKNKCLsp12887268", "WSSKNKCLsp12140271",
                "WSSKNKCLsp12887270", "WSSKNKCLsp12140273",
            ],
        },
    ],
    "Stomach": [
        {
            "files": [
                "Itai_GSE203612/raw_h5ad/GSM6177607_NYU_GIST1_Vis_raw.h5ad",
                "Itai_GSE203612/raw_h5ad/GSM6177609_NYU_GIST2_Vis_raw.h5ad",
            ]
        },
        {"glob": "Gastric_GSE251950/raw_h5ad/"},
    ],
    "Thyroid": [
        {"glob": "Thyroid_SNU/raw_h5ad/"},
    ],
}

# "Major lineage membership (for *_enriched proportion columns)" -----------
LINEAGE_SUBTYPES = {
    "Epi": [
        "Epi_Cycle", "Epi_Luminal", "Epi_Hormonal", "Epi_Xenobiotic", "Epi_TNF",
        "Epi_Interferon", "Epi_Glandular", "Epi_Squamous", "Epi_Hypox-Stress",
        "Epi_cEMT", "Epi_ER-Stress", "Epi_Basal", "Epi_MHCII", "Epi_Hypox-Adapt",
        "Epi_Ciliated", "Epi_Oxphos-Ion", "Epi_Oxphos-Metal",
    ],
    "T": [
        "CD8_Trm", "CD8_Tem/Trm", "CD16-_NK", "CD8_Tex",
        "CD4_Tem/Effector", "CD4_Treg", "CD8_Tem/Temra", "ILC3",
        "CD4_Tfh", "CD4_Tcm/Tn", "CD4_CXCL13+_T", "MAIT",
        "γδ_T", "CD4_Th17", "NK", "CD8_Tcm/Tn",
        "CD4_Unassigned_T", "CD8_Cycling_Tex", "CD16+_NK", "CD8_ISG15+_T", "CD4_Cycling_Treg",
        "CD4_Th1",
    ],
    "B": [
        "NR4A2+_Bn", "Bn", "Bmem", "Bgc", "Plasma",
        "Cycling_Bgc", "ABC", "Plasmablast",
    ],
    "Mye": [
        "CXCL3+_Macro", "CD14+_Mono", "CD1C+_cDC", "FOLR2+_Macro", "HSPA6+_Macro",
        "SPP1+_cycMacro", "CD16+_Mono", "Macro", "LILRA4+_pDC", "CLEC9A+_cDC",
        "TNFSF10+_Macro", "SPP1+_CXCL3+_Macro", "SPP1+_Macro", "SPP1+_MT1+_Macro",
        "LAMP3+_cDC", "Mast",
    ],
    "EC": [
        "Venous_EC", "Capillary_EC", "Arterial_EC", "Fibrosis_PGF+_Tip_EC",
        "Venous_iEC", "PGF+_Tip_EC", "Lymphatic_EC", "HMGB2+_cycEC", "Tip_EC",
        "SPRY1+_EC", "CCL2+_iEC", "HMOX1+_iEC", "CD14+_EC", "TMEM100+_EC",
    ],
    "Fib": [
        "CXCL14+_mCAF", "PI16+_iCAF", "iCAF", "Fibroblast", "apCAF", "vCAF",
        "IL6+_iCAF", "pnCAF", "mCAF", "ISG15+_mCAF", "tCAF", "HSP+_tCAF",
        "cyc_mCAF", "myoCAF",
    ],
    "Mu": ["Pericyte", "SMC"],
}

MODEL_PARAMS = dict(N_cells_per_location=30, detection_alpha=20)
TRAIN_PARAMS = dict(max_epochs=30000, train_size=1, use_gpu=True)


# "Data assembly" -----------------------------------------------------------
def build_data_list(organ):
    data_list = []
    for entry in ORGAN_DATASETS[organ]:
        if "files" in entry:
            data_list.extend(entry["files"])
            continue
        path = entry["glob"]
        files = sorted(os.listdir(os.path.join(INPUT_DIR, path)))
        files = [path + f for f in files]
        if "keep_identifiers" in entry:
            ids = set(entry["keep_identifiers"])
            files = [f for f in files if any(i in f for i in ids)]
        data_list.extend(files)
    return data_list


def drop_low_quality_slides(data_list):
    # Drop a slide if QC filtering removes >50% of spots, or if it has <15000 genes.
    keep = []
    for rel_path in data_list:
        adata_tmp = sc.read(os.path.join(INPUT_DIR, rel_path))
        raw_ncell = adata_tmp.shape[0]

        sc.pp.filter_cells(adata_tmp, min_counts=500)
        adata_tmp = adata_tmp[adata_tmp.obs["pct_counts_mt"] < 30]
        filter_ncell = adata_tmp.shape[0]

        bad_quality = (raw_ncell - filter_ncell) / raw_ncell > 0.5
        bad_feature = adata_tmp.shape[1] < 15000
        if bad_quality or bad_feature:
            reason = "quality" if bad_quality else "feature"
            print(f"  dropped ({reason}): {rel_path}")
            continue
        keep.append(rel_path)
    return keep


def concat_visium_slides(data_list):
    adata_vis = None
    for rel_path in data_list:
        adata_tmp = sc.read(os.path.join(INPUT_DIR, rel_path))
        sample_geo = adata_tmp.obs["Sample_GEO"].unique()[0]
        adata_tmp.obs_names = [f"{sample_geo}@{n}" for n in adata_tmp.obs_names]

        sc.pp.filter_cells(adata_tmp, min_counts=500)
        adata_tmp = adata_tmp[adata_tmp.obs["pct_counts_mt"] < 30]

        if adata_vis is None:
            adata_vis = adata_tmp
            continue
        adata_vis = anndata.concat(
            [adata_vis, adata_tmp], join="inner", uns_merge="unique", index_unique=None
        )
    return adata_vis


# "Deconvolution" -------------------------------------------------------------
def run_deconvolution(adata_vis, gpu):
    inf_aver = pd.read_csv(f"{REF_RUN_NAME}/inf_aver.csv", index_col=0)

    adata_vis.obs_names_make_unique()
    adata_vis.var["SYMBOL"] = adata_vis.var_names

    mt_gene = [gene.startswith("MT-") for gene in adata_vis.var["SYMBOL"]]
    adata_vis.var["MT_gene"] = mt_gene
    adata_vis.obsm["MT"] = adata_vis[:, adata_vis.var["MT_gene"].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var["MT_gene"].values]

    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample_id")
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver, **MODEL_PARAMS
    )

    torch.cuda.set_device(gpu)
    mod.train(batch_size=TRAIN_BATCH_SIZE, **TRAIN_PARAMS)

    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={"num_samples": 1000, "batch_size": mod.adata.n_obs, "use_gpu": True},
    )
    return adata_vis, mod


def save_qc_plots(mod, run_name2):
    os.makedirs(run_name2, exist_ok=True)

    mod.plot_history(1000)
    plt.legend(labels=["full data training"])
    plt.savefig(f"{run_name2}/ELBO_loss.png", bbox_inches="tight")
    plt.close()

    mod.plot_QC()
    plt.savefig(f"{run_name2}/QC.png", bbox_inches="tight")
    plt.close()

    mod.plot_spatial_QC_across_batches()
    plt.savefig(f"{run_name2}/spatial_QC_across_batches.png", bbox_inches="tight")
    plt.close()


# "Post-processing" -----------------------------------------------------------
def add_race_info(adata_vis):
    info = pd.read_csv(SLIDE_INFO_PATH, sep="\t")
    data_geo = adata_vis.obs["Data_GEO"].unique()

    country = [
        info.loc[info["Data_GEO"] == g, "Country"].iloc[0]
        if not info[info["Data_GEO"] == g].empty else None
        for g in data_geo
    ]
    race = [
        info.loc[info["Data_GEO"] == g, "Race"].iloc[0]
        if not info[info["Data_GEO"] == g].empty else None
        for g in data_geo
    ]
    race_info = pd.DataFrame({"Data_GEO": data_geo, "Country": country, "Race": race})

    merged_obs = pd.merge(adata_vis.obs, race_info, on="Data_GEO", how="left")
    merged_obs.index = adata_vis.obs_names
    adata_vis.obs = merged_obs
    return adata_vis


def compute_lineage_proportions(adata_vis):
    factor_names = np.array(
        [x.replace(" ", "_") for x in adata_vis.uns["mod"]["factor_names"]]
    )
    adata_vis.uns["mod"]["factor_names"] = factor_names
    adata_vis.obs[factor_names] = adata_vis.obsm["q05_cell_abundance_w_sf"]

    row_sum = adata_vis.obs[factor_names].sum(axis=1)
    for cell_type in factor_names:
        adata_vis.obs[cell_type] = adata_vis.obs[cell_type] / row_sum

    for lineage, subtypes in LINEAGE_SUBTYPES.items():
        adata_vis.obs[f"{lineage}_enriched"] = adata_vis.obs[subtypes].sum(axis=1)

    return adata_vis


# "Main" ------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--organ", choices=sorted(ORGAN_CANCER_TYPE), help="Organ to run")
    parser.add_argument("--gpu", type=int, default=DEFAULT_GPU, help="CUDA device index")
    parser.add_argument("--list-organs", action="store_true", help="List available organs and exit")
    args = parser.parse_args()

    if args.list_organs:
        for organ, cancer_type in sorted(ORGAN_CANCER_TYPE.items()):
            print(f"{organ}\t{cancer_type}")
        return
    if args.organ is None:
        parser.error("--organ is required (or use --list-organs)")

    organ = args.organ
    gpu = args.gpu

    # matches the fixed value used in every notebook (not tied to torch.cuda.set_device)
    os.environ["THEANO_FLAGS"] = "device=cuda0,floatX=float32,force_device=True"
    os.chdir(PROJECT_DIR)

    run_name2 = f"{RUN_NAME}cell2location_map_{organ}_30000epoch"
    print(f"## Organ: {organ} ({ORGAN_CANCER_TYPE[organ]}), GPU {gpu}")

    print("# Building Visium dataset list")
    data_list = build_data_list(organ)
    print(f"  {len(data_list)} candidate slides")

    print("# QC filtering slides")
    data_list = drop_low_quality_slides(data_list)
    print(f"  {len(data_list)} slides retained")

    print("# Concatenating slides")
    adata_vis = concat_visium_slides(data_list)
    print(f"  {adata_vis.shape}")

    print("# Running cell2location deconvolution")
    adata_vis, mod = run_deconvolution(adata_vis, gpu)

    print("# Saving model and QC plots")
    mod.save(run_name2, overwrite=True)
    save_qc_plots(mod, run_name2)

    adata_vis.obs.columns = adata_vis.obs.columns.astype(str)
    adata_vis.var.columns = adata_vis.var.columns.astype(str)
    adata_vis.write(f"{run_name2}/sp.h5ad")

    print("# Adding race/country info")
    adata_vis = add_race_info(adata_vis)

    print("# Computing lineage proportions")
    adata_vis = compute_lineage_proportions(adata_vis)

    print("# Saving final outputs")
    adata_vis.obs.to_csv(f"{run_name2}/Obs.csv")
    adata_vis.to_df().to_csv(f"{run_name2}/Counts.csv")
    adata_vis.write(f"{run_name2}/sp_race.h5ad")

    print(f"## Done: {run_name2}")


if __name__ == "__main__":
    main()
