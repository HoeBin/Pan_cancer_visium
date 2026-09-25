# 분석 목적
#   - Figure 5a, 5b를 재현한다 (pan-cancer 재현성, LR z-score).
#
# 분석 흐름
#   1. Fig.5a 암종별 주요 cell-type pair z-score heatmap과 재현 암종 수 barplot을 그린다
#      (Source_Data_Fig5.xlsx 사용).
#   2. Fig.5b T cell-B cell 및 Endothelial-Mural LR pair z-score heatmap을 그린다
#      (Source_Data_Fig5.xlsx의 Fig.5b top 시트 사용, 아래 검증 참고).
#
# 주요 출력
#   - Fig5a.pdf, Fig5b.pdf
#
# 출력 위치
#   - Output/Figure5/
#
# 참고: Fig.5c(UECA/HNCA/BRCA/OVCA/LUCA 샘플의 spatial co-enrichment 이미지)는
# Figure5/Fig5c_spatial_lr_overlap.py(외부 Visium h5ad 필요)가, Fig.5d(scRNA-seq/Visium
# LR strength 암종별 density/상관관계)는 Figure5/Fig5d_sc_vs_visium_scatter.py가 담당한다.
#
# 참고(Fig5b 데이터 출처): 예전 Source_Data_Fig5.xlsx의 "Fig.5b top/bottom" 시트는 실제
# 출판된 Fig5b와 다른 LR pair 목록을 담고 있어(FST-BMPR2 등 논문 유전자 없음) 원본 분석과
# 대조한 별도 CSV(Verified_Source_Data/)를 썼다. 시트가 수정되어 그 CSV와 값이 일치하므로
# (988행, LR pair 64개, z-score 차이 < 1e-15) CSV를 제거하고 시트를 직접 읽는다. 시트가 다시
# 바뀌어 행 수·LR pair 수·암종 수가 아래 기대값과 다르면 조용히 다른 그림을 그리지 않도록
# 오류로 멈춘다. 의도적으로 갱신한 경우 fig5b_expected를 새 값으로 고친다.

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

### Source Data 검증 ####
# group_pair별 행 수 / LR pair 수 (Fig.5b: T cell-B cell 37개, Endothelial-Mural 39개 LR pair x 13개 암종)
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
  lr_order <- unique(x$lr_pair)          # Figure의 행 순서 (파일에 이미 그 순서로 기록됨)
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
