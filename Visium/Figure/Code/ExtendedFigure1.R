# Purpose
#   - Reproduces Extended Fig.1c and 1f from Source_Data_Extended_Fig1.xlsx
#     (MSigDB Hallmark 2020 pathway enrichment per Myeloid/Fibroblast subtype).
#
# Workflow
#   1. Arrange the Enrichr results (Term x cluster_name) into a -log10(Adjusted P-value) matrix.
#   2. Draw a heatmap that also marks significance levels (*, **, ***) in the cells.
#
# Main outputs
#   - ExtFig1c.pdf (Myeloid), ExtFig1f.pdf (Fibroblast)
#
# Output location
#   - Output/ExtendedFigure1/
#
# Note: Extended Fig.1a-b and 1d-e (scRNA-seq UMAP, marker gene dot plot) are embedding/dot plots without
# Source Data, so they are excluded.

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ComplexHeatmap); library(circlize); library(viridisLite)
})
script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(script_path), "common/config.R"))

out_dir  <- file.path(output_dir, "ExtendedFigure1/Tables")
plot_dir <- file.path(output_dir, "ExtendedFigure1/Plots")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

sig_symbol <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}

plot_hallmark <- function(sheet, panel, title) {
  dat <- read_source("Source_Data_Extended_Fig1.xlsx", sheet)
  dat <- dat[!is.na(dat$Term), ]
  names(dat)[names(dat) == "-log10(Adjusted P-value)'"] <- "neglog10padj"
  write.csv(dat, file.path(out_dir, paste0(panel, "_hallmark_enrichment.csv")), row.names = FALSE)

  mat <- xtabs(neglog10padj ~ Term + cluster_name, data = dat)
  mat[mat == 0] <- 0
  star_dat <- dat %>% select(Term, cluster_name, `Adjusted P-value`) %>%
    mutate(symbol = sig_symbol(`Adjusted P-value`))
  star_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  for (i in seq_len(nrow(star_dat))) {
    star_mat[star_dat$Term[i], star_dat$cluster_name[i]] <- star_dat$symbol[i]
  }

  # Saturate the color scale at 5 (-log10(padj)), same as the original legend.
  mat_capped <- pmin(mat, 5)
  col_fun <- colorRamp2(c(0, 2.5, 5), magma(3))
  ht <- Heatmap(mat_capped, name = "-log10(padj)", col = col_fun,
    column_title = title,
    cluster_rows = TRUE, cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9), column_names_rot = 45,
    heatmap_legend_param = list(at = c(0, 1.25, 2.5, 3.75, 5)),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(star_mat[i, j], x, y, gp = gpar(fontsize = 8, col = "black"))
    })
  pdf(file.path(plot_dir, paste0(panel, ".pdf")), width = 8, height = 8)
  draw(ht)
  dev.off()
}

plot_hallmark("ED Fig.1c", "ExtFig1c", "Myeloid")
plot_hallmark("ED Fig.1f", "ExtFig1f", "Fibroblast")

message("Extended Figure 1 complete: ", dirname(out_dir))
