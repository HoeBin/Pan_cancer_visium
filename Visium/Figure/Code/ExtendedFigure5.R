# 분석 목적
#   - Source_Data_Extended_Fig5.xlsx로부터 Extended Fig.5a를 재현한다
#     (major cell-type pair co-localization의 compartment별 비교).
#
# 분석 흐름
#   1. Fig.5a 21개 major cell-type pair의 compartment별(Malignant/Boundary/Normal)
#      median Jaccard score 분포를 grouped boxplot과 pair 표시 패널로 그린다.
#
# 주요 출력
#   - ExtFig5a.pdf
#
# 출력 위치
#   - Output/ExtendedFigure5/
#
# 참고: Ext.Fig.5b(Epi-Epi/Epi-TME/TME-TME 세 그룹의 compartment 간 비교)는
# Code/Figure3/05_ExtDataFig5b_7d_pair_categories.py가 담당한다 (같은
# Source_Data_Extended_Fig5.xlsx의 "ED Fig.5b" 시트). 해당 스크립트는 같은 슬라이드가
# 세 compartment 모두에 매칭되는 데이터임을 반영해 paired Wilcoxon(+슬라이드 매칭)을
# 쓰며, 이 R 스크립트가 이전에 쓰던 `stat_compare_means(..., comparisons=...)`는
# paired=TRUE 지정이 없어 unpaired Wilcoxon으로 계산되던 차이가 있었다.

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})
script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(script_path), "common/config.R"))

out_dir  <- file.path(output_dir, "ExtendedFigure5/Tables")
plot_dir <- file.path(output_dir, "ExtendedFigure5/Plots")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# "Extended Fig 5a: major cell-type pair by compartment" -----------------------
a5 <- read_source("Source_Data_Extended_Fig5.xlsx", "ED Fig.5a")
a5 <- a5[!is.na(a5$Feat_pair), ]

# 21개 pair는 값으로 정렬하지 않는다. legend 순서(Epithelial -> T cell -> B cell ->
# Myeloid -> Endothelial -> Fibroblast -> Mural)의 상삼각(i<j) 조합을 그대로 나열하는
# 고정된 순서다 (원본 그림에서 각 anchor 블록 내부도 celltype_order를 그대로 따름).
all_pairs <- unique(a5$Feat_pair)
pair_c1 <- sub("-.*$", "", all_pairs)
pair_c2 <- sub("^.*?-", "", all_pairs)
anchor_labels <- celltype_labels[celltype_order]
pair_order <- character(0)
for (i in seq_along(anchor_labels)) {
  for (j in seq_along(anchor_labels)) {
    if (j <= i) next
    match_idx <- which((pair_c1 == anchor_labels[i] & pair_c2 == anchor_labels[j]) |
                        (pair_c1 == anchor_labels[j] & pair_c2 == anchor_labels[i]))
    pair_order <- c(pair_order, all_pairs[match_idx])
  }
}
a5$Feat_pair <- factor(a5$Feat_pair, levels = pair_order)
a5$Region <- factor(a5$Region, levels = mal_order)
write.csv(a5, file.path(out_dir, "ExtFig5a_plot_data.csv"), row.names = FALSE)

p5a_top <- ggplot(a5, aes(Feat_pair, median_J_comp, fill = Region)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = mal_colors, name = "Malignancy") +
  base_theme +
  labs(x = NULL, y = "Median Jaccard score") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Pan_cancer_visium/Visium/Figure/Supplementary_Figure15.R(p15a_matrix)의 원본
## 패턴을 따른다: y축은 "top/bottom" 같은 임의 값이 아니라 실제 7개 cell type이며,
## 각 pair는 자신을 이루는 두 cell type의 고정된 행(celltype_order 순서, 위 Epithelial
## ~ 아래 Mural)에 점을 찍고 grey 선으로 잇는 진짜 UpSet-matrix 구조다.
pair_dots <- bind_rows(
  a5 %>% distinct(Feat_pair) %>% mutate(celltype = sub("-.*$", "", Feat_pair)),
  a5 %>% distinct(Feat_pair) %>% mutate(celltype = sub("^.*-", "", Feat_pair))
)
pair_dots$celltype <- factor(pair_dots$celltype, levels = rev(anchor_labels))
p5a_bottom <- ggplot(pair_dots, aes(x = Feat_pair, y = celltype, group = Feat_pair)) +
  geom_line(color = "grey40", linewidth = 1.2) +
  geom_point(aes(color = celltype), size = 3) +
  scale_color_manual(values = celltype_colors_spaced, name = "Cell Type", drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = "Cell-type pair") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_text(angle = 90, size = 12),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey85"))

pdf(file.path(plot_dir, "ExtFig5a.pdf"), width = 13, height = 6.5)
print(p5a_top / p5a_bottom + plot_layout(heights = c(5, 1), guides = "collect"))
dev.off()

message("Extended Figure 5a complete: ", dirname(out_dir))
