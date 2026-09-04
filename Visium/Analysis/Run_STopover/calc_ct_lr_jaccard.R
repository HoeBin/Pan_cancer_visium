# Purpose
#   - Combine STopover's ligand-receptor (LR) pair spatial co-localization results
#     (jaccard_composite_lr.csv, output of run_stopover.py) with cell2location-based
#     cell type spatial connected components (Comb_CC_*) to compute the Jaccard
#     index between "where a given cell type is present" and "where a given LR
#     pair is active".
#   - Consolidates the 5 calculation blocks that were repeated in
#     STopver_result/CT_LR_Jaccard.R for each scope (Global major lineage / Sub
#     fine-grained subtype) and compartment (all spots / Mal / Bdy / Normal)
#     combination into a single script.
#
# Analysis flow
#   1. Load per-slide STopover LR results and the celltype+LR merged obs table
#      (the Sub scope additionally merges spatial_final_celltype_obs.csv into
#      spatial_final_mod_obs.csv).
#   2. Binarize spots to O/X based on whether the Comb_CC_* value is non-zero for
#      the cell type pair (group1, group2) and the LR pair (lig, rec), then compute
#      the Jaccard index from the resulting contingency table.
#   3. Save results across the Global (7 major lineages, 21 pairs) / Sub-Interest
#      (Myeloid subtype vs. T/Endothelial/Fibroblast subtype) / Sub-All (28 major
#      lineage pair groups expanded to subtype level) scopes, each combined with
#      an all-spots / compartment (Mal/Bdy/Normal) filter.
#
# Main outputs
#   - GroupPair_LR_Jaccard_each_sample*.csv per scope/cell-type-pair/compartment
#     combination
#
# Output location
#   - /path/to/output/lr_jaccard/global_all_ct/
#   - /path/to/output/lr_jaccard/sub_all_ct/
#
# Note: the paths below are example directory names. Edit them for your
# environment before running. lr_dir must point to the same location as the
# "global" scope's new_output_root in run_stopover.py (jaccard_composite_lr.csv).

# Libraries and paths ------------------------------------------------------------
library(stringr)
library(parallel)

lr_dir <- "/path/to/output/stopover_lr/global_all/"
global_ct_dir <- "/path/to/input/stopover_celltype/global_all"
sub_ct_dir <- "/path/to/input/stopover_celltype/sub_all"
global_out_dir <- "/path/to/output/lr_jaccard/global_all_ct/"
sub_out_dir <- "/path/to/output/lr_jaccard/sub_all_ct/"
dir.create(global_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sub_out_dir, recursive = TRUE, showWarnings = FALSE)

# Cell type subtype membership (major lineage -> fine subtype) -------------------
major_lineages <- c("Epithelial", "T cell", "B cell", "Myeloid", "Endothelial", "Fibroblast", "Mural")

s_Epithelial <- c(
  "Epi_Cycle", "Epi_Luminal", "Epi_Hormonal", "Epi_Xenobiotic", "Epi_TNF",
  "Epi_Interferon", "Epi_Glandular", "Epi_Squamous", "Epi_Hypox-Stress",
  "Epi_cEMT", "Epi_ER-Stress", "Epi_Basal", "Epi_MHCII", "Epi_Hypox-Adapt",
  "Epi_Ciliated", "Epi_Oxphos-Ion", "Epi_Oxphos-Metal"
)
s_Tcell <- c(
  "CD8_Trm", "CD8_Tem/Trm", "CD16-_NK", "CD8_Tex",
  "CD4_Tem/Effector", "CD4_Treg", "CD8_Tem/Temra", "ILC3",
  "CD4_Tfh", "CD4_Tcm/Tn", "CD4_CXCL13+_T", "MAIT",
  "γδ_T", "CD4_Th17", "NK", "CD8_Tcm/Tn",
  "CD4_Unassigned_T", "CD8_Cycling_Tex", "CD16+_NK", "CD8_ISG15+_T", "CD4_Cycling_Treg",
  "CD4_Th1"
)
s_Bcell <- c(
  "NR4A2+_Bn", "Bn", "Bmem", "Bgc", "Plasma",
  "Cycling_Bgc", "ABC", "Plasmablast"
)
s_Myeloid <- c(
  "CXCL3+_Macro", "CD14+_Mono", "CD1C+_cDC", "FOLR2+_Macro", "HSPA6+_Macro",
  "SPP1+_cycMacro", "CD16+_Mono", "Macro", "LILRA4+_pDC", "CLEC9A+_cDC",
  "TNFSF10+_Macro", "SPP1+_CXCL3+_Macro", "SPP1+_Macro", "SPP1+_MT1+_Macro",
  "LAMP3+_cDC", "Mast"
)
s_Endothelial <- c(
  "Venous_EC", "Capillary_EC", "Arterial_EC", "Fibrosis_PGF+_Tip_EC",
  "Venous_iEC", "PGF+_Tip_EC", "Lymphatic_EC", "HMGB2+_cycEC", "Tip_EC",
  "SPRY1+_EC", "CCL2+_iEC", "HMOX1+_iEC", "CD14+_EC", "TMEM100+_EC"
)
s_Fibroblast <- c(
  "CXCL14+_mCAF", "PI16+_iCAF", "iCAF", "Normal_Fibroblast", "apCAF", "vCAF",
  "IL6+_iCAF", "pnCAF", "mCAF", "ISG15+_mCAF", "tCAF", "HSP+_tCAF",
  "cyc_mCAF", "myoCAF"
)
s_Mural <- c("SMC", "Pericyte")

lineage_subtypes <- list(
  s_Epithelial = s_Epithelial, s_Tcell = s_Tcell, s_Bcell = s_Bcell, s_Myeloid = s_Myeloid,
  s_Fibroblast = s_Fibroblast, s_Endothelial = s_Endothelial, s_Mural = s_Mural
)

# Core Jaccard calculation --------------------------------------------------------
### Jaccard index of a contingency table ####
jaccard_index <- function(mat) {
  jaccard_matrix <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = list(rownames(mat), colnames(mat)))
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      intersection <- mat[i, j]
      union <- sum(mat[i, ]) + sum(mat[, j]) - intersection
      jaccard_matrix[i, j] <- if (union > 0) intersection / union else NA
    }
  }
  jaccard_matrix
}

### Cell type pair x LR pair Jaccard for one slide (optionally within one compartment) ####
compute_slide_jaccard <- function(obs, lr, combinations, group_sep, loc = NA) {
  if (!is.na(loc)) obs <- obs[obs$LocationNew == loc, ]

  result_list <- list()
  idx <- 1
  for (j in 1:length(combinations)) {
    celltype1 <- combinations[[j]][1]
    celltype2 <- combinations[[j]][2]

    for (k in 1:nrow(lr)) {
      a <- lr[k, ]
      obs_tmp <- obs[c(
        "cell_id", "LocationNew",
        paste0("Comb_CC_", celltype1), paste0("Comb_CC_", celltype2),
        paste0("Comb_CC_", a$Feat_1), paste0("Comb_CC_", a$Feat_2)
      )]
      colnames(obs_tmp) <- c("cellid", "LocationNew", "group1", "group2", "lig", "rec")
      obs_tmp$group1 <- ifelse(obs_tmp$group1 != 0, "O", "X")
      obs_tmp$group2 <- ifelse(obs_tmp$group2 != 0, "O", "X")
      obs_tmp$lig <- ifelse(obs_tmp$lig != 0, "O", "X")
      obs_tmp$rec <- ifelse(obs_tmp$rec != 0, "O", "X")
      obs_tmp$group_pair <- paste0(obs_tmp$group1, "-", obs_tmp$group2)
      obs_tmp$lr_pair <- paste0(obs_tmp$lig, "-", obs_tmp$rec)

      contingency_table <- table(obs_tmp$group_pair, obs_tmp$lr_pair)
      jaccard_matrix <- jaccard_index(contingency_table)
      jaccard_value <- if ("O-O" %in% rownames(jaccard_matrix) && "O-O" %in% colnames(jaccard_matrix)) {
        jaccard_matrix["O-O", "O-O"]
      } else {
        0
      }

      tmp_df <- data.frame(
        group_pair = paste0(celltype1, group_sep, celltype2),
        lr_pair = paste0(a$Feat_1, "-", a$Feat_2),
        cancer_type = as.character(unique(obs$cancer_type)),
        sample_id = as.character(unique(obs$sample_id)),
        jaccard = jaccard_value
      )
      if (!is.na(loc)) tmp_df$LocationNew <- loc

      result_list[[idx]] <- tmp_df
      idx <- idx + 1
    }
  }
  do.call(rbind, result_list)
}

### Run the calculation over every slide in a folder, in parallel ####
run_ct_lr_jaccard <- function(folder_list, ct_dir, combinations, group_sep, merge_celltype_obs, loc = NA, num_cores = 20) {
  process_folder <- function(i) {
    file <- folder_list[i]
    obs <- read.csv(paste0(ct_dir, "/", file, "/spatial_final_mod_obs.csv"), row.names = 1, check.names = FALSE)
    obs$cell_id <- rownames(obs)

    if (merge_celltype_obs) {
      obs_ct <- read.csv(paste0(ct_dir, "/", file, "/spatial_final_celltype_obs.csv"), row.names = 1, check.names = FALSE)
      obs_ct$cell_id <- rownames(obs_ct)
      ct_cols <- colnames(obs_ct)[str_detect(colnames(obs_ct), "Comb_CC_")]
      obs <- cbind(obs, obs_ct[, ct_cols])
    }

    lr <- read.csv(paste0(lr_dir, "/", file, "/jaccard_composite_lr.csv"), check.names = FALSE)
    compute_slide_jaccard(obs, lr, combinations, group_sep, loc)
  }
  results <- mclapply(1:length(folder_list), process_folder, mc.cores = num_cores)
  do.call(rbind, results)
}

# Global scope: 7 major lineages, 21 pairs ----------------------------------------
folder_list_global <- list.files(global_ct_dir)
global_combinations <- combn(major_lineages, 2, simplify = FALSE)

### All spots ####
df_total <- run_ct_lr_jaccard(folder_list_global, global_ct_dir, global_combinations,
  group_sep = "-", merge_celltype_obs = FALSE, num_cores = 40
)
write.csv(df_total, paste0(global_out_dir, "GroupPair_LR_Jaccard_each_sample.csv"), quote = FALSE, row.names = FALSE)

### Per compartment (Mal/Bdy/Normal) ####
for (loc in c("Mal", "Bdy", "Normal")) {
  df_total <- run_ct_lr_jaccard(folder_list_global, global_ct_dir, global_combinations,
    group_sep = "-", merge_celltype_obs = FALSE, loc = loc, num_cores = 40
  )
  write.csv(df_total, paste0(global_out_dir, "GroupPair_LR_Jaccard_each_sample_", loc, ".csv"), quote = FALSE, row.names = FALSE)
}

# Sub scope: fine subtype pairs ----------------------------------------------------
folder_list_sub <- list.files(sub_ct_dir)
sub_all_combinations <- combn(c(s_Epithelial, s_Tcell, s_Bcell, s_Myeloid, s_Fibroblast, s_Endothelial, s_Mural), 2, simplify = FALSE)

### Interest pairs: Myeloid subtype vs Tcell/Endothelial/Fibroblast subtype, all spots ####
interest_combinations <- Filter(function(p) {
  (p[1] %in% s_Myeloid && p[2] %in% c(s_Tcell, s_Endothelial, s_Fibroblast)) ||
    (p[2] %in% s_Myeloid && p[1] %in% c(s_Tcell, s_Endothelial, s_Fibroblast))
}, sub_all_combinations)

df_total <- run_ct_lr_jaccard(folder_list_sub, sub_ct_dir, interest_combinations,
  group_sep = "@", merge_celltype_obs = TRUE, num_cores = 20
)
write.csv(df_total, paste0(sub_out_dir, "GroupPair_LR_Jaccard_each_sample_Interest.csv"), quote = FALSE, row.names = FALSE)

### All major-lineage pair groups (incl. self pairs), exploded to subtype level ####
lineage_names <- names(lineage_subtypes)
lineage_pair_groups <- c(combn(lineage_names, 2, simplify = FALSE), lapply(lineage_names, function(x) c(x, x)))

for (pair_group in lineage_pair_groups) {
  lineage1 <- pair_group[1]
  lineage2 <- pair_group[2]
  label1 <- gsub("s_", "", lineage1)
  label2 <- gsub("s_", "", lineage2)

  combinations <- Filter(function(p) {
    (p[1] %in% lineage_subtypes[[lineage1]] && p[2] %in% lineage_subtypes[[lineage2]]) ||
      (p[2] %in% lineage_subtypes[[lineage2]] && p[1] %in% lineage_subtypes[[lineage1]])
  }, sub_all_combinations)

  #### All spots ####
  df_total <- run_ct_lr_jaccard(folder_list_sub, sub_ct_dir, combinations,
    group_sep = "@", merge_celltype_obs = TRUE, num_cores = 20
  )
  write.csv(df_total, paste0(sub_out_dir, "GroupPair_LR_Jaccard_each_sample_", label1, "@", label2, ".csv"), quote = FALSE, row.names = FALSE)

  #### Per compartment (Mal/Bdy/Normal) ####
  for (loc in c("Mal", "Bdy", "Normal")) {
    df_total <- run_ct_lr_jaccard(folder_list_sub, sub_ct_dir, combinations,
      group_sep = "@", merge_celltype_obs = TRUE, loc = loc, num_cores = 20
    )
    write.csv(df_total, paste0(sub_out_dir, "GroupPair_LR_Jaccard_each_sample_", label1, "@", label2, "_", loc, ".csv"), quote = FALSE, row.names = FALSE)
  }
}
