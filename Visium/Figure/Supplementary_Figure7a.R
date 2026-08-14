# 분석 목적
#   - 세포 유형 비율–marker 발현 상관관계 중간 결과에서 Supplementary Fig.7a를 재현한다.
#
# 분석 흐름
#   1. 사전 계산된 correlation RDS를 불러온다.
#   2. 7개 대표 marker와 global cell type을 선택한다.
#   3. violin/boxplot과 marker annotation을 생성한다.
#   4. plot data 및 marker별 요약 결과를 CSV로 저장한다.
#
# 주요 출력
#   - Supplementary_Fig7a.pdf
#   - marker별 correlation 원자료 및 요약 CSV
#
# 출력 위치
#   - Figure/Output/Supplementary_Figure7/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
root_dir <- dirname(script_path)
source(file.path(root_dir, "common/config.R"))
out_dir <- file.path(root_dir, "Output/Supplementary_Figure7/Tables")
plot_dir <- file.path(root_dir, "Output/Supplementary_Figure7/Plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- file.path(PROJECT_DIR, "Cell2Location_SubNum4/CorPropCounts/Output/Large_CorPropCounts.rds")
gene_list <- c("CDH1", "CD3E", "CD79A", "CD68", "RAMP2", "COL1A1", "ACTA2")
cell_levels <- c("Epithelial", "Tcell", "Bcell", "Myeloid", "Endothelial", "Fibroblast", "Mural")
marker_map <- data.frame(gene = gene_list, type = cell_levels, stringsAsFactors = FALSE)
colors <- c(Epithelial = "red", Tcell = "blue", Bcell = "green", Myeloid = "darkorchid", Endothelial = "magenta", Fibroblast = "orange", Mural = "#00EEEE")

# "Load and filter intermediate correlations" --------------------------------
result <- readRDS(input_rds)
result <- as.data.frame(result)
result <- result[result$cancer_type != "Esophagus", ]
type_map <- c(Epi = "Epithelial", T = "Tcell", B = "Bcell", Mye = "Myeloid", EC = "Endothelial", Fib = "Fibroblast", Mu = "Mural")
result$type <- unname(type_map[result$type])
result <- result[!is.na(result$cor) & result$gene %in% gene_list & result$type %in% cell_levels, ]
result$gene <- factor(result$gene, levels = gene_list)
result$type <- factor(result$type, levels = cell_levels)

write.csv(result, file.path(out_dir, "Supplementary_Fig7a_plot_data_nonESCA.csv"), row.names = FALSE)
write.csv(result, file.path(out_dir, "Supplementary_Fig7a_plot_data.csv"), row.names = FALSE)
summary_table <- result %>%
  group_by(type, gene) %>%
  summarise(n = n(), median_correlation = median(cor, na.rm = TRUE), mean_correlation = mean(cor, na.rm = TRUE), fraction_p_lt_0_05 = mean(p_value < 0.05, na.rm = TRUE), .groups = "drop")
write.csv(summary_table, file.path(out_dir, "Supplementary_Fig7a_key_results_nonESCA.csv"), row.names = FALSE)
write.csv(summary_table, file.path(out_dir, "Supplementary_Fig7a_key_results.csv"), row.names = FALSE)

# "Reproduce plot" ------------------------------------------------------------
p_main <- ggplot(result, aes(gene, cor, fill = gene)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  geom_boxplot(width = 0.2, outlier.size = 1, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = unname(colors[marker_map$type])) +
  scale_y_continuous(expand = c(0, 0), position = "right", breaks = c(-0.5, 0, 0.5)) +
  facet_grid(rows = vars(type), scales = "free", switch = "y") +
  theme_bw() + labs(x = NULL, y = "Correlation", title = "Correlation between cell type proportion\nand marker gene expression") +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"), plot.title = element_text(hjust = 0.5), strip.background = element_blank(), strip.text = element_text(size = 20), strip.text.y.left = element_text(angle = 0), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 18))

marker_map$gene <- factor(marker_map$gene, levels = gene_list)
p_annotation <- ggplot(marker_map, aes(gene, 1, fill = type)) +
  geom_tile() + scale_fill_manual(values = colors) + scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(direction = "vertical", title = NULL, keyheight = 0.5, nrow = 1)) +
  theme_void() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1, size = 20, color = "black"))

pdf(file.path(plot_dir, "Supplementary_Fig7a.pdf"), width = 9, height = 12)
print(plot_grid(p_main, p_annotation, ncol = 1, rel_heights = c(0.85, 0.15), align = "v", axis = "lr"))
dev.off()

message("Supplementary Figure 7a complete: ", dirname(out_dir))
