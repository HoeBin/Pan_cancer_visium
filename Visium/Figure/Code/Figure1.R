# Purpose
#   - Reproduces the per-cancer-type numbers of scRNA-seq/Visium samples, cells and spots from Source_Data_Fig1.xlsx (Fig.1a).
#
# Workflow
#   1. Read the Fig.1a sheet and tabulate the sample size per cancer type.
#   2. Draw bar plots for scRNA-seq (samples/cells) and Visium (samples/spots) separately.
#   3. Combine the four metrics into one patchwork figure and save it.
#
# Main outputs
#   - Fig1a.pdf (bar plots of scRNA-seq/Visium samples, cells and spots)
#   - Fig1a_sample_summary.csv
#
# Output location
#   - Output/Figure1/
#
# Note: the original Fig.1a is a circular (fan) bar plot built around a BioRender human silhouette, and
# Fig.1b is entirely a BioRender workflow schematic. The copyrighted BioRender artwork and the manual
# layout are not reproduced in code; only the numerical data in Fig.1a is reproduced here as bar plots
# with the same information content. Fig.1b is excluded because it has no data to reproduce.

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(scales)
})
script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(script_path), "common/config.R"))

out_dir  <- file.path(output_dir, "Figure1/Tables")
plot_dir <- file.path(output_dir, "Figure1/Plots")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# "Figure 1a" -------------------------------------------------------------
fig1a <- read_source("Source_Data_Fig1.xlsx", "Fig.1a")
fig1a <- fig1a[fig1a$Cancer_type != "Total", ]
fig1a$Cancer_type <- factor(fig1a$Cancer_type, levels = cancer_order[cancer_order %in% fig1a$Cancer_type])
write.csv(fig1a, file.path(out_dir, "Fig1a_sample_summary.csv"), row.names = FALSE)

bar_panel <- function(dat, y, ylab, fill) {
  ggplot(dat, aes(x = Cancer_type, y = .data[[y]])) +
    geom_col(fill = fill, color = "black", width = 0.7) +
    scale_y_continuous(labels = comma) +
    labs(x = NULL, y = ylab) +
    base_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

p_sc_samples  <- bar_panel(fig1a, "sc_n_samples", "scRNA-seq samples", "#8B0000")
p_sc_cells    <- bar_panel(fig1a, "sc_n_cells",   "scRNA-seq cells",   "#CD5C5C")
p_vis_samples <- bar_panel(fig1a, "vis_n_samples", "Visium samples",   "#00008B")
p_vis_spots   <- bar_panel(fig1a, "vis_n_spots",   "Visium spots",     "#6495ED")

p1a <- (p_sc_samples | p_sc_cells) / (p_vis_samples | p_vis_spots) +
  plot_annotation(title = "Figure 1a: sample, cell and spot counts by cancer type")

pdf(file.path(plot_dir, "Fig1a.pdf"), width = 11, height = 8)
print(p1a)
dev.off()

message("Figure 1 complete: ", dirname(out_dir))
