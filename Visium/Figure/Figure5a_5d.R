# 분석 목적
#   - Figure 5a와 5d를 plot-ready 중간 CSV에서 재현한다.
#
# 분석 흐름
#   1. Fig5a cell-pair recurrence/z-score CSV를 읽어 heatmap-bar plot을 생성한다.
#   2. Fig5d scRNA-seq–Visium LR strength CSV를 읽어 암종별 density plot을 생성한다.
#   3. plot 입력 및 암종별 주요 통계를 CSV로 저장한다.
#
# 주요 출력
#   - Fig5a 및 Fig5d 재현 PDF
#   - cell-pair recurrence와 scRNA–Visium concordance 결과 CSV
#
# 출력 위치
#   - Figure/Output/Figure5/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
root_dir <- dirname(script_path)
source(file.path(root_dir, "common/config.R"))
source_dir <- file.path(PROJECT_DIR, "Cell2Location_SubNum4/Fig5_code/Code")
out_dir  <- file.path(root_dir, "Output/Figure5/Tables")
plot_dir <- file.path(root_dir, "Output/Figure5/Plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# "Figure 5a" -----------------------------------------------------------------
old_wd <- getwd()
setwd(source_dir)
source("Fig5a_plot_from_python_output_v2.R", local = new.env(parent = globalenv()))
setwd(old_wd)

fig5a_input <- fread(file.path(source_dir, "figures/Fig5a_R_input/top5_cellpair_heatmap_bar_input_for_R_zscore_by_cancer_type.csv"))
fig5a_summary <- fread(file.path(source_dir, "figures/Fig5a_R_input/top5_cellpair_pair_summary.csv"))
fig5a_input <- fig5a_input[cancer_type != "ESCA"]
fig5a_input[, mask_flag := tolower(as.character(mask)) %in% c("true", "t", "1")]
fig5a_summary <- fig5a_input[, .(
  n_cancer_topk = sum(mask_flag),
  mean_topk_j = mean(raw_J_comp_agg[mask_flag], na.rm = TRUE),
  mean_rank = mean(row_order[mask_flag], na.rm = TRUE)
), by = cell_pair]
fig5a_input[, mask_flag := NULL]
write.csv(fig5a_input, file.path(out_dir, "Fig5a_plot_data_nonESCA.csv"), row.names = FALSE)
write.csv(fig5a_summary, file.path(out_dir, "Fig5a_key_results_cellpair_recurrence_nonESCA.csv"), row.names = FALSE)
write.csv(fig5a_input, file.path(out_dir, "Fig5a_plot_data.csv"), row.names = FALSE)
write.csv(fig5a_summary, file.path(out_dir, "Fig5a_key_results_cellpair_recurrence.csv"), row.names = FALSE)

fig5a_pdf <- file.path(source_dir, "figures/Fig5a_R_output/top5_cellpair_recurrence_bar_heatmap_aligned_zscore_by_cancer_type_ComplexHeatmap_R_nonESCA.pdf")
file.copy(fig5a_pdf, file.path(plot_dir, "Fig5a.pdf"), overwrite = TRUE)

# "Figure 5d" -----------------------------------------------------------------
setwd(source_dir)
source("Fig5d_plot_from_python_output.R", local = new.env(parent = globalenv()))
setwd(old_wd)

fig5d_input <- fread(file.path(source_dir, "figures/sc_vs_visium_bothsig_plot_input_for_R_all_target_pairs.csv"))
fig5d_input <- fig5d_input[cancer_type != "ESCA"]
write.csv(fig5d_input, file.path(out_dir, "Fig5d_plot_data_nonESCA.csv"), row.names = FALSE)

# "Figure 5b: LR z-score heatmaps" -------------------------------------------
make_lr_zscore <- function(group_name, top_n = 40) {
  x <- fig5d_input[group_pair == group_name, .(
    strength = median(sc_strength, na.rm = TRUE)
  ), by = .(lr_pair, cancer_type)]
  x[, group_pair := group_name]
  x[, z_score := if (.N > 1 && stats::sd(strength, na.rm = TRUE) > 0) as.numeric(scale(strength)) else 0,
    by = lr_pair]
  selected <- x[, .(max_abs_z = max(abs(z_score), na.rm = TRUE)), by = lr_pair][
    order(-max_abs_z), head(lr_pair, top_n)
  ]
  x[lr_pair %in% selected]
}
fig5b <- rbindlist(list(
  make_lr_zscore("T cell-B cell"),
  make_lr_zscore("Tcell-Bcell"),
  make_lr_zscore("B cell-T cell"),
  make_lr_zscore("Endothelial-Mural")
), use.names = TRUE, fill = TRUE)
fig5b <- unique(fig5b)
write.csv(fig5b, file.path(out_dir, "Fig5b_key_results_LR_zscores_nonESCA.csv"), row.names = FALSE)
write.csv(fig5b, file.path(out_dir, "Fig5b_key_results_LR_zscores.csv"), row.names = FALSE)

if (nrow(fig5b) > 0) {
  suppressPackageStartupMessages({library(ComplexHeatmap); library(circlize)})
  pdf(file.path(plot_dir, "Fig5b.pdf"), width = 8, height = 9)
  for (group_name in unique(fig5b$group_pair)) {
    x <- as.data.frame(fig5b[group_pair == group_name, .(lr_pair, cancer_type, z_score)])
    mat <- xtabs(z_score ~ lr_pair + cancer_type, data = x)
    draw(Heatmap(mat, name = "z-score", column_title = group_name,
      col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_names_gp = grid::gpar(fontsize = 6), column_names_rot = 45))
  }
  dev.off()
}

safe_cor <- function(x, y, method) {
  keep <- complete.cases(x, y)
  if (sum(keep) < 3 || length(unique(x[keep])) < 2 || length(unique(y[keep])) < 2) return(NA_real_)
  cor(x[keep], y[keep], method = method)
}
fig5d_summary <- fig5d_input[, .(
  n_lr_pairs = .N,
  n_group_pairs = uniqueN(group_pair),
  median_sc_strength = median(sc_strength, na.rm = TRUE),
  median_visium_strength = median(vis_strength, na.rm = TRUE),
  pearson_r = safe_cor(sc_strength, vis_strength, "pearson"),
  spearman_rho = safe_cor(sc_strength, vis_strength, "spearman")
), by = cancer_type]
write.csv(fig5d_summary, file.path(out_dir, "Fig5d_key_results_by_cancer.csv"), row.names = FALSE)

fig5d_pdf <- file.path(source_dir, "figures/sc_vs_visium_by_target_pair_R/sc_vs_visium_bothsig_scatter_by_cancer_with_contour_paircolor_labeled_R_all_target_pairs_nonESCA.pdf")
file.copy(fig5d_pdf, file.path(plot_dir, "Fig5d.pdf"), overwrite = TRUE)

message("Figure 5a and 5d complete: ", dirname(out_dir))
