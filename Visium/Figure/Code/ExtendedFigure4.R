# Purpose
#   - Reproduces Extended Fig.4a-d from Source_Data_Extended_Fig4.xlsx
#     (CNV score/spot count per compartment, correlation of the CNV score with cell-type proportions).
#
# Workflow
#   1. Draw median-CNV and spot-count boxplots with slide-matched paired Wilcoxon tests (Fig.4a-b).
#   2. Draw a bar plot of the correlation between global cell-type proportions and the CNV score (Fig.4c).
#   3. Draw a bar plot of the correlation between subtype proportions and the CNV score, colored by lineage (Fig.4d).
#
# Main outputs
#   - ExtFig4a.pdf, ExtFig4b.pdf, ExtFig4c.pdf, ExtFig4d.pdf
#   - ExtFig4ab_paired_statistics.csv (sample-matched two-sided P values)
#
# Output location
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
paired_statistics <- function(dat, y, panel) {
  # Match by slide ID explicitly; Excel row order must not determine pairing.
  if (anyNA(dat$sample_id) || anyDuplicated(dat[c("sample_id", "Malignancy")]))
    stop(panel, ": missing sample IDs or duplicate slide-compartment rows")
  wide <- dat %>% select(sample_id, Malignancy, all_of(y)) %>%
    pivot_wider(names_from = Malignancy, values_from = all_of(y))
  if (!all(mal_order %in% names(wide))) stop(panel, ": missing compartments")
  comparisons <- list(c("Malignant", "Boundary"), c("Boundary", "Normal"),
                      c("Malignant", "Normal"))
  bind_rows(lapply(comparisons, function(pair) {
    keep <- is.finite(wide[[pair[1]]]) & is.finite(wide[[pair[2]]])
    if (sum(keep) < 2L) stop(panel, ": fewer than two matched slides")
    test <- wilcox.test(wide[[pair[1]]][keep], wide[[pair[2]]][keep],
                        paired = TRUE, alternative = "two.sided", exact = FALSE, correct = TRUE)
    data.frame(panel = panel, variable = y, group1 = pair[1], group2 = pair[2],
               n_pairs = sum(keep), test = "wilcoxon_signed_rank", p_value = test$p.value)
  }))
}

plot_compartment <- function(dat, y, title, ylab, stats) {
  dat$Malignancy <- factor(dat$Malignancy, levels = mal_order)
  yrange <- range(dat[[y]], finite = TRUE)
  span <- max(diff(yrange), 1)
  stats$y.position <- yrange[2] + span * c(0.08, 0.18, 0.28)
  stats$label <- sprintf("p = %.2e", stats$p_value)
  ggplot(dat, aes(Malignancy, .data[[y]])) +
    geom_jitter(width = 0.15, size = 0.6, color = "grey60", alpha = 0.4) +
    geom_boxplot(aes(color = Malignancy), fill = NA, outlier.shape = NA, linewidth = 0.8) +
    stat_pvalue_manual(stats, label = "label", tip.length = 0.01) +
    scale_color_manual(values = mal_colors, guide = "none") +
    base_theme +
    labs(title = title, x = NULL, y = ylab)
}

a4 <- read_source("Source_Data_Extended_Fig4.xlsx", "ED Fig.4a")
a4 <- a4[!is.na(a4$Malignancy), ]
write.csv(a4, file.path(out_dir, "ExtFig4a_cnv_score.csv"), row.names = FALSE)
a4_stats <- paired_statistics(a4, "median_cnv", "ExtFig4a")
pdf(file.path(plot_dir, "ExtFig4a.pdf"), width = 4.5, height = 5)
print(plot_compartment(a4, "median_cnv", "CNV score", "CNV", a4_stats))
dev.off()

b4 <- read_source("Source_Data_Extended_Fig4.xlsx", "ED Fig.4b")
b4 <- b4[!is.na(b4$Malignancy), ]
write.csv(b4, file.path(out_dir, "ExtFig4b_spot_count.csv"), row.names = FALSE)
b4_stats <- paired_statistics(b4, "n", "ExtFig4b")
write.csv(bind_rows(a4_stats, b4_stats), file.path(out_dir, "ExtFig4ab_paired_statistics.csv"), row.names = FALSE)
pdf(file.path(plot_dir, "ExtFig4b.pdf"), width = 4.5, height = 5)
print(plot_compartment(b4, "n", "Number of spot", "Count", b4_stats))
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
# Use the cairo_pdf device because subtype names contain Greek letters (gamma-delta T).
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
