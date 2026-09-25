# 분석 목적
#   - Source_Data_Extended_Fig4.xlsx로부터 Extended Fig.4a-d를 재현한다
#     (compartment별 CNV score/spot 수, CNV score와 cell-type proportion의 상관관계).
#
# 분석 흐름
#   1. Fig.4a-b compartment별 CNV score, spot 수 boxplot을 그린다.
#   2. Fig.4c global cell-type proportion과 CNV score의 상관관계 bar plot을 그린다.
#   3. Fig.4d subtype proportion과 CNV score의 상관관계 bar plot을 lineage 색상으로 그린다.
#
# 주요 출력
#   - ExtFig4a.pdf, ExtFig4b.pdf, ExtFig4c.pdf, ExtFig4d.pdf
#
# 출력 위치
#   - Output/ExtendedFigure4/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggpubr)
})
script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(script_path), "common/config.R"))

out_dir  <- file.path(output_dir, "ExtendedFigure4/Tables")
plot_dir <- file.path(output_dir, "ExtendedFigure4/Plots")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# "Extended Fig 4a-b: CNV score / spot count by compartment" ------------------
plot_compartment <- function(dat, y, title, ylab) {
  dat$Malignancy <- factor(dat$Malignancy, levels = mal_order)
  ggplot(dat, aes(Malignancy, .data[[y]])) +
    geom_jitter(width = 0.15, size = 0.6, color = "grey60", alpha = 0.4) +
    geom_boxplot(aes(color = Malignancy), fill = NA, outlier.shape = NA, linewidth = 0.8) +
    stat_compare_means(method = "wilcox.test", comparisons = mal_comparisons, label = "p.format") +
    scale_color_manual(values = mal_colors, guide = "none") +
    base_theme +
    labs(title = title, x = NULL, y = ylab)
}

a4 <- read_source("Source_Data_Extended_Fig4.xlsx", "ED Fig.4a")
a4 <- a4[!is.na(a4$Malignancy), ]
write.csv(a4, file.path(out_dir, "ExtFig4a_cnv_score.csv"), row.names = FALSE)
pdf(file.path(plot_dir, "ExtFig4a.pdf"), width = 4.5, height = 5)
print(plot_compartment(a4, "mean_cnv", "CNV score", "CNV"))
dev.off()

b4 <- read_source("Source_Data_Extended_Fig4.xlsx", "ED Fig.4b")
b4 <- b4[!is.na(b4$Malignancy), ]
write.csv(b4, file.path(out_dir, "ExtFig4b_spot_count.csv"), row.names = FALSE)
pdf(file.path(plot_dir, "ExtFig4b.pdf"), width = 4.5, height = 5)
print(plot_compartment(b4, "n", "Number of spot", "Count"))
dev.off()

# "Extended Fig 4c: global proportion vs CNV score correlation" ---------------
c4 <- read_source("Source_Data_Extended_Fig4.xlsx", "ED Fig.4c")
c4 <- c4[!is.na(c4$celltype), ] %>% arrange(desc(correlation))
c4$celltype <- factor(c4$celltype, levels = c4$celltype)
write.csv(c4, file.path(out_dir, "ExtFig4c_global_correlation.csv"), row.names = FALSE)
p4c <- ggplot(c4, aes(celltype, correlation, fill = celltype)) +
  geom_col(color = "black") +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  base_theme +
  labs(title = "Correlation with global proportion and CNV score", x = NULL, y = "Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file.path(plot_dir, "ExtFig4c.pdf"), width = 5, height = 5)
print(p4c)
dev.off()

# "Extended Fig 4d: sub proportion vs CNV score correlation" ------------------
subtype_lineage <- load_subtype_lineage()
d4 <- read_source("Source_Data_Extended_Fig4.xlsx", "ED Fig.4d")
d4 <- d4[!is.na(d4$celltype), ] %>% arrange(desc(correlation))
d4$celltype <- factor(d4$celltype, levels = d4$celltype)
d4$lineage <- subtype_lineage[as.character(d4$celltype)]
write.csv(d4, file.path(out_dir, "ExtFig4d_sub_correlation.csv"), row.names = FALSE)
# subtype 이름에 그리스 문자(gamma-delta T)가 있어 cairo_pdf 장치를 사용한다.
p4d <- ggplot(d4, aes(celltype, correlation, fill = lineage)) +
  geom_col(color = NA) +
  scale_fill_manual(values = celltype_colors_spaced, name = "Lineage") +
  base_theme +
  labs(title = "Correlation with sub proportion and CNV score", x = NULL, y = "Correlation") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
cairo_pdf(file.path(plot_dir, "ExtFig4d.pdf"), width = 16, height = 6)
print(p4d)
dev.off()

message("Extended Figure 4 complete: ", dirname(out_dir))
