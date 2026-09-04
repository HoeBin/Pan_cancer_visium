# 분석 목적
#   - Figure 5a, 5b를 재현한다 (pan-cancer 재현성, LR z-score).
#
# 분석 흐름
#   1. Fig.5a 암종별 주요 cell-type pair z-score heatmap과 재현 암종 수 barplot을 그린다
#      (Source_Data_Fig5.xlsx 사용).
#   2. Fig.5b T cell-B cell 및 Endothelial-Mural LR pair z-score heatmap을 그린다
#      (Verified_Source_Data/ 사용, 아래 참고).
#
# 주요 출력
#   - Fig5a.pdf, Fig5b.pdf
#
# 출력 위치
#   - Output/Figure5/
#
# 참고: Fig.5c(UECA/HNCA/BRCA/OVCA/LUCA 샘플의 spatial co-enrichment 이미지)와 Fig.5d
# (scRNA-seq/Visium LR strength 암종별 density/상관관계)는 이번 범위에서 제외한다.
#
# 참고(Fig5b 데이터 출처): 공식 배포된 Source_Data_Fig5.xlsx의 "Fig.5b top/bottom"
# 시트는 실제 출판된 Fig5b의 LR pair 목록과 전혀 다른 값(예: FST-BMPR2, THBS2-CD36 등
# 실제 논문 유전자가 하나도 없음)이 담겨 있어 사용할 수 없다. raw per-sample STopover
# Jaccard(GroupPair_LR_Jaccard_each_sample.csv)와 원본 분석 스크립트(CT_LR_Jaccard.R,
# celltype_order 21개 pair 중 "T cell-B cell"=13p, "Endothelial-Mural"=39p)를 대조해
# 이 대안이 실제 published Fig5b와 유전자 목록·순서·값이 정확히 일치함을 확인했다
# (Code_data_figure/R/Figure5b_LR_SpecificCommon.R). 그 결과를
# Verified_Source_Data/Fig5b_key_results_LR_zscores.csv로 로컬에 포함해 이 저장소가
# 여전히 self-contained하게 실행되도록 했다.

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(ComplexHeatmap); library(circlize)
})
script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(script_path), "common/config.R"))

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
fig5b_verified <- read.csv(file.path(dirname(script_path), "Verified_Source_Data/Fig5b_key_results_LR_zscores.csv"),
                            stringsAsFactors = FALSE)
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

message("Figure 5 complete: ", dirname(out_dir))
