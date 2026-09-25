# Purpose
#   - Reproduces Extended Fig.5a from Source_Data_Extended_Fig5.xlsx
#     (compartment-wise comparison of major cell-type pair co-localization).
#
# Workflow
#   1. Draw the distribution of the median Jaccard score of the 21 major cell-type pairs (Fig.5a) per compartment
#      (Malignant/Boundary/Normal) as a grouped boxplot with a pair-indicator panel.
#
# Main outputs
#   - ExtFig5a.pdf
#
# Output location
#   - Output/ExtendedFigure5/
#
# Note: Ext.Fig.5b (comparison of the three groups Epi-Epi/Epi-TME/TME-TME across compartments) is handled
# by Code/Figure3/ExtDataFig5b_7d_pair_categories.py (the "ED Fig.5b" sheet of the same
# Source_Data_Extended_Fig5.xlsx). That script uses a paired Wilcoxon test (with slide matching), reflecting
# that the same slide is matched in all three compartments, whereas the `stat_compare_means(..., comparisons=...)`
# previously used in this R script did not specify paired=TRUE and therefore computed an unpaired Wilcoxon test.

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

# The 21 pairs are not sorted by value. They are a fixed order that simply lists the upper-triangle (i<j)
# combinations of the legend order (Epithelial -> T cell -> B cell -> Myeloid -> Endothelial -> Fibroblast -> Mural)
# (in the original figure, the order inside each anchor block also follows celltype_order).
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

## Follows the original pattern of Pan_cancer_visium/Visium/Figure/Supplementary_Figure15.R (p15a_matrix): the y-axis is not
## an arbitrary value such as "top/bottom" but the actual 7 cell types, and each pair is a true UpSet-matrix
## structure: points are placed on the fixed rows of the two cell types that make it up (celltype_order,
## Epithelial at the top to Mural at the bottom) and joined by a grey line.
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
