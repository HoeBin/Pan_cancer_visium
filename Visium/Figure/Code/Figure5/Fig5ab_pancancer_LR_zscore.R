# Purpose
#   - Reproduces Figures 5a and 5b (pan-cancer reproducibility, LR z-score).
#
# Workflow
#   1. Draw the per-cancer-type z-score heatmap of the top cell-type pairs and the barplot of the number of
#      reproducing cancer types for Fig.5a (uses Source_Data_Fig5.xlsx).
#   2. Draw the z-score heatmaps of the T cell-B cell and Endothelial-Mural LR pairs for Fig.5b
#      (uses the "Fig.5b top" sheet of Source_Data_Fig5.xlsx; see the validation below).
#
# Main outputs
#   - Fig5a.pdf, Fig5b.pdf
#
# Output location
#   - Output/Figure5/
#
# Note: Fig.5c (spatial co-enrichment images of UECA/HNCA/BRCA/OVCA/LUCA samples) is handled by
# Figure5/Fig5c_spatial_lr_overlap.py (needs an external Visium h5ad), and Fig.5d (per-cancer-type density/
# correlation of scRNA-seq/Visium LR strength) by Figure5/Fig5d_sc_vs_visium_scatter.py.
#
# Note (source of the Fig5b data): the "Fig.5b top/bottom" sheets of an earlier Source_Data_Fig5.xlsx contained an
# LR pair list different from the published Fig5b (none of the paper's genes such as FST-BMPR2), so a separate
# CSV (Verified_Source_Data/), checked against the original analysis, was used instead. The sheet has since been
# corrected and agrees with that CSV (988 rows, 64 LR pairs, z-score difference < 1e-15), so the CSV was removed and
# the sheet is read directly. If the sheet changes again and the number of rows, LR pairs or cancer types differs
# from the expected values below, the script stops with an error instead of silently drawing a different figure.
# If the change is intentional, update fig5b_expected to the new values.

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(ComplexHeatmap); library(circlize)
})
script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir    <- dirname(dirname(script_path))
source(file.path(code_dir, "common/config.R"))

out_dir  <- file.path(output_dir, "Figure5/Tables")
plot_dir <- file.path(output_dir, "Figure5/Plots")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# "Figure 5a: pan-cancer cell-pair recurrence" ---------------------------------
a5 <- read_source("Source_Data_Fig5.xlsx", "Fig.5a")
a5 <- a5[!is.na(a5$cell_pair), ]
write.csv(a5, file.path(out_dir, "Fig5a_plot_data.csv"), row.names = FALSE)

pair_row <- a5 %>% distinct(cell_pair, row_order, n_cancer_topk) %>% arrange(row_order)
cancer_col <- a5 %>% distinct(cancer_type, col_order) %>% arrange(col_order)
mat5a <- xtabs(value ~ cell_pair + cancer_type, data = a5)
mat5a <- mat5a[pair_row$cell_pair, cancer_col$cancer_type]

fill_min <- unique(a5$fill_min)[1]; fill_mid <- unique(a5$fill_mid)[1]; fill_max <- unique(a5$fill_max)[1]
col_fun5a <- colorRamp2(c(fill_min, fill_mid, fill_max), c("#2166AC", "white", "#B2182B"))
bar_anno <- rowAnnotation(
  "n cancer types" = anno_barplot(pair_row$n_cancer_topk, gp = gpar(fill = "grey40"), width = unit(2.5, "cm")))
ht5a <- Heatmap(mat5a, name = "z-score", col = col_fun5a,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 11), column_names_gp = gpar(fontsize = 11), column_names_rot = 45,
  right_annotation = bar_anno,
  heatmap_legend_param = list(title = unique(a5$fill_label)[1]))
pdf(file.path(plot_dir, "Fig5a.pdf"), width = 9, height = 6)
draw(ht5a)
dev.off()

# "Figure 5b: LR pair z-score heatmaps" ----------------------------------------
fig5b_verified <- read_source("Source_Data_Fig5.xlsx", "Fig.5b top")
fig5b_verified <- fig5b_verified[!is.na(fig5b_verified$lr_pair), c("lr_pair", "cancer_type", "strength", "group_pair", "z_score")]

### Source Data validation ####
# Rows / LR pairs per group_pair (Fig.5b: 37 T cell-B cell and 39 Endothelial-Mural LR pairs x 13 cancer types)
fig5b_expected <- data.frame(
  group_pair = c("T cell-B cell", "Endothelial-Mural"),
  n_rows     = c(481L, 507L),
  n_lr_pair  = c(37L, 39L),
  n_cancer   = c(13L, 13L)
)
fig5b_expected_total <- list(n_rows = 988L, n_lr_pair = 64L)
fig5b_problems <- character()
for (i in seq_len(nrow(fig5b_expected))) {
  x <- fig5b_verified[fig5b_verified$group_pair == fig5b_expected$group_pair[i], ]
  got <- c(n_rows = nrow(x), n_lr_pair = length(unique(x$lr_pair)), n_cancer = length(unique(x$cancer_type)))
  want <- unlist(fig5b_expected[i, c("n_rows", "n_lr_pair", "n_cancer")])
  if (!all(got == want)) {
    fig5b_problems <- c(fig5b_problems, sprintf("%s: rows/LR pair/cancer = %s (expected %s)",
                                                fig5b_expected$group_pair[i], paste(got, collapse = "/"), paste(want, collapse = "/")))
  }
}
if (nrow(fig5b_verified) != fig5b_expected_total$n_rows)
  fig5b_problems <- c(fig5b_problems, sprintf("total rows = %d (expected %d)", nrow(fig5b_verified), fig5b_expected_total$n_rows))
if (length(unique(fig5b_verified$lr_pair)) != fig5b_expected_total$n_lr_pair)
  fig5b_problems <- c(fig5b_problems, sprintf("total LR pairs = %d (expected %d)", length(unique(fig5b_verified$lr_pair)), fig5b_expected_total$n_lr_pair))
if (!all(fig5b_verified$cancer_type %in% cancer_order))
  fig5b_problems <- c(fig5b_problems, paste("unknown cancer_type:", paste(setdiff(fig5b_verified$cancer_type, cancer_order), collapse = ", ")))
if (any(duplicated(fig5b_verified[, c("lr_pair", "cancer_type", "group_pair")])))
  fig5b_problems <- c(fig5b_problems, "duplicated (lr_pair, cancer_type, group_pair) rows")
if (anyNA(fig5b_verified$z_score))
  fig5b_problems <- c(fig5b_problems, "NA in z_score")
if (length(fig5b_problems) > 0)
  stop("Source_Data_Fig5.xlsx 'Fig.5b top' does not match the expected published Fig.5b:\n  ",
       paste(fig5b_problems, collapse = "\n  "), call. = FALSE)
write.csv(fig5b_verified, file.path(out_dir, "Fig5b_plot_data.csv"), row.names = FALSE)

make_lr_heatmap <- function(group_label) {
  x <- fig5b_verified[fig5b_verified$group_pair == group_label, ]
  lr_order <- unique(x$lr_pair)          # row order of the figure (already recorded in that order in the file)
  col_order <- intersect(cancer_order, unique(x$cancer_type))
  mat <- xtabs(z_score ~ lr_pair + cancer_type, data = x)
  mat <- mat[lr_order, col_order]
  Heatmap(mat, name = "z-score", column_title = group_label,
    col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 10), column_names_rot = 45)
}
pdf(file.path(plot_dir, "Fig5b.pdf"), width = 8, height = 9)
draw(make_lr_heatmap("T cell-B cell"))
draw(make_lr_heatmap("Endothelial-Mural"))
dev.off()

message("Figure 5a-b complete: ", dirname(out_dir))
