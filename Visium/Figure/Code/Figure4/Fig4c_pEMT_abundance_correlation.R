# Purpose
#   - Reproduces Figure 4c (pEMT-abundance correlation rank plot) from Source_Data_Fig4.xlsx.
#
# Workflow
#   1. Draw the rank of the Spearman correlation of each cell subtype with the pEMT score as a scatter plot.
#
# Main outputs
#   - Fig4c.pdf
#
# Output location
#   - Output/Figure4/
#
# Note: Fig.4a/b (signature scores per compartment) are handled by Figure4/Fig4ab_signature_panels.py
# (recomputed with paired Wilcoxon + Friedman, 268 slide-matched). Fig.4d (LR pair rank),
# Fig.4f (hazard ratio forest), Fig.4g (Kaplan-Meier) and Fig.4e (PACA sample spatial multi-panel) are handled
# by the Python scripts in Figure4/.

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
})
script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir    <- dirname(dirname(script_path))
source(file.path(code_dir, "common/config.R"))

out_dir  <- file.path(output_dir, "Figure4/Tables")
plot_dir <- file.path(output_dir, "Figure4/Plots")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

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

message("Figure 4c complete: ", dirname(out_dir))
