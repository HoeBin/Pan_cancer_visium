# 분석 목적
#   - Source_Data_Fig4.xlsx로부터 Figure 4a-c를 재현한다 (compartment별 gene-set score).
#
# 분석 흐름
#   1. Fig.4a-b compartment별 signature score(Epithelial/Mesenchymal/pEMT) boxplot을 그린다.
#   2. Fig.4c pEMT-abundance correlation rank plot을 그린다.
#
# 주요 출력
#   - Fig4a_Epithelial.pdf, Fig4a_Mesenchymal.pdf, Fig4b.pdf, Fig4c.pdf
#
# 출력 위치
#   - Output/Figure4/
#
# 참고: Fig.4d(LR pair rank), Fig.4f(hazard ratio forest), Fig.4g(Kaplan-Meier)와
# Fig.4e(PACA 샘플 spatial multi-panel)는 이번 범위에서 제외한다.

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggpubr); library(ggrepel)
})
script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(script_path), "common/config.R"))

out_dir  <- file.path(output_dir, "Figure4/Tables")
plot_dir <- file.path(output_dir, "Figure4/Plots")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# "Figure 4a-b: compartment signature scores" ----------------------------------
plot_score <- function(dat, title) {
  dat$Malignancy <- factor(dat$Malignancy, levels = mal_order)
  ggplot(dat, aes(Malignancy, value)) +
    geom_jitter(width = 0.15, size = 0.6, color = "grey60", alpha = 0.4) +
    geom_boxplot(aes(color = Malignancy), fill = NA, outlier.shape = NA, linewidth = 0.8) +
    stat_compare_means(method = "wilcox.test", comparisons = mal_comparisons, label = "p.format") +
    scale_color_manual(values = mal_colors, guide = "none") +
    base_theme +
    labs(title = title, x = NULL, y = "Median score")
}

a4 <- read_source("Source_Data_Fig4.xlsx", "Fig.4a")
a4 <- a4[!is.na(a4$signature), ]
for (sig in unique(a4$signature)) {
  x <- a4[a4$signature == sig, ]
  x_long <- pivot_longer(x, c("Malignant", "Boundary", "Normal"), names_to = "Malignancy", values_to = "value")
  write.csv(x_long, file.path(out_dir, paste0("Fig4a_", sig, "_plot_data.csv")), row.names = FALSE)
  pdf(file.path(plot_dir, paste0("Fig4a_", sig, ".pdf")), width = 4, height = 5)
  print(plot_score(x_long, sig))
  dev.off()
}

b4 <- read_source("Source_Data_Fig4.xlsx", "Fig.4b")
b4 <- b4[!is.na(b4$signature), ]
b4_long <- pivot_longer(b4, c("Malignant", "Boundary", "Normal"), names_to = "Malignancy", values_to = "value")
write.csv(b4_long, file.path(out_dir, "Fig4b_pEMT_plot_data.csv"), row.names = FALSE)
pdf(file.path(plot_dir, "Fig4b.pdf"), width = 4, height = 5)
print(plot_score(b4_long, "pEMT"))
dev.off()

# "Figure 4c: pEMT-abundance correlation" --------------------------------------
c4 <- read_source("Source_Data_Fig4.xlsx", "Fig.4c")
c4 <- c4[!is.na(c4$celltype), ]
c4$global_celltype <- recode(c4$global_celltype, `T cell` = "Tcell", `B cell` = "Bcell")
c4$global_celltype <- factor(c4$global_celltype, levels = celltype_order)
write.csv(c4, file.path(out_dir, "Fig4c_pEMT_abundance_correlations.csv"), row.names = FALSE)
top_label <- bind_rows(slice_head(c4, n = 6), slice_tail(c4, n = 10))
p4c <- ggplot(c4, aes(rank, spearman_rho)) +
  geom_point(aes(color = global_celltype), size = 2) +
  geom_text_repel(data = top_label, aes(label = celltype), size = 3.5, max.overlaps = 20) +
  scale_color_manual(values = celltype_colors, labels = celltype_labels, name = "Cell type") +
  base_theme +
  labs(title = "pEMT-abundance correlation", x = "Rank", y = "Spearman rho") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
pdf(file.path(plot_dir, "Fig4c.pdf"), width = 8, height = 5)
print(p4c)
dev.off()

message("Figure 4 complete: ", dirname(out_dir))
