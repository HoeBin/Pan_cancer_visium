# Purpose
#   - Reproduces the abundance-correlation heatmaps (with dot annotation) between major cell types for Fig.3a (all 13
#     cancer types) and Extended Data Fig.7a (BRCA downsampled to 30 slides).
#
# Workflow
#   1. Read the input table (celltype_1, celltype_2, correlation).
#      - Fig.3a   : the Fig.3a sheet of Source_Data_Fig3.xlsx (Paper_info_260928) by default;
#                   --fig3a-source=csv falls back to the frozen Input/fig3a_source_data.csv
#      - Ext.Fig.7a: the ED Fig.7a sheet of Source_Data_Extended_Fig7.xlsx
#   2. Convert the long table into a 7x7 correlation matrix and check its symmetry and value range.
#   3. Draw a ComplexHeatmap whose rows/columns are marked by cell-type color dots.
#
# Main outputs
#   - fig3a_global_correlation.pdf, extfig7a_global_correlation.pdf
#   - source_data_fig3a.csv, source_data_extfig7a.csv
#
# Output location
#   - Output/Figure3/, Output/ExtendedFigure7/
#
# Note: Source_Data_Fig3.xlsx's Fig.3a sheet was corrected in the Paper_info_260928 snapshot (2026-09-28) and now
# matches the frozen Input/fig3a_source_data.csv (abundance-based correlation, ESCA excluded, 13 cancer types) to
# within floating-point precision, so the script now reads the xlsx sheet by default. The script always prints the
# max absolute difference between the two so a future Source Data update that silently drifts is caught; pass
# --fig3a-source=csv to force the frozen snapshot instead of the live sheet.

# Libraries and paths ---------------------------------------------------------
suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(grid)
})
script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir    <- dirname(dirname(script_path))
source(file.path(code_dir, "common/config.R"))

fig3a_source <- "xlsx"
for (arg in commandArgs(trailingOnly = TRUE)) {
  if (grepl("^--fig3a-source=", arg)) fig3a_source <- sub("^--fig3a-source=", "", arg)
}
if (!fig3a_source %in% c("csv", "xlsx")) stop("--fig3a-source must be csv or xlsx")
message("Fig.3a source: ", fig3a_source)

fig3a_csv <- file.path(code_dir, "Figure3/Input/fig3a_source_data.csv")

# Output locations per panel (out_dir: tables, plot_dir: plots)
panels <- list(
  fig3a    = list(out_dir = file.path(output_dir, "Figure3/Tables"),
                  plot_dir = file.path(output_dir, "Figure3/Plots")),
  extfig7a = list(out_dir = file.path(output_dir, "ExtendedFigure7/Tables"),
                  plot_dir = file.path(output_dir, "ExtendedFigure7/Plots"))
)
for (panel in names(panels)) {
  dir.create(panels[[panel]]$out_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(panels[[panel]]$plot_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load source data ------------------------------------------------------------
src <- list()
fig3a_csv_data  <- read.csv(fig3a_csv, stringsAsFactors = FALSE)
fig3a_xlsx_data <- read_source("Source_Data_Fig3.xlsx", "Fig.3a")
fig3a_xlsx_data <- fig3a_xlsx_data[!is.na(fig3a_xlsx_data$celltype_1), c("celltype_1", "celltype_2", "correlation")]
# Always compare the frozen snapshot against the live sheet so a future Source Data update that silently
# drifts from the frozen snapshot is caught, regardless of which one is actually plotted below.
chk <- merge(fig3a_csv_data, fig3a_xlsx_data, by = c("celltype_1", "celltype_2"), suffixes = c("_csv", "_xlsx"))
message("Fig.3a: frozen Input/fig3a_source_data.csv vs Source_Data_Fig3.xlsx max |diff| = ",
        signif(max(abs(chk$correlation_csv - chk$correlation_xlsx)), 3),
        " (", nrow(chk), " / ", nrow(fig3a_csv_data), " pairs matched)")
src$fig3a <- if (fig3a_source == "xlsx") fig3a_xlsx_data else fig3a_csv_data
src$extfig7a <- read_source("Source_Data_Extended_Fig7.xlsx", "ED Fig.7a")

# Draw heatmaps ---------------------------------------------------------------
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

for (panel in names(panels)) {
  df <- src[[panel]]
  df <- df[!is.na(df$celltype_1), c("celltype_1", "celltype_2", "correlation")]
  write.csv(df, file.path(panels[[panel]]$out_dir, paste0("source_data_", panel, ".csv")), row.names = FALSE)

  ### Long table -> correlation matrix ####
  mat <- matrix(NA_real_, nrow = length(celltype_order), ncol = length(celltype_order),
                dimnames = list(celltype_order, celltype_order))
  for (i in seq_len(nrow(df))) {
    mat[df$celltype_1[i], df$celltype_2[i]] <- df$correlation[i]
  }
  if (anyNA(mat)) stop(panel, ": correlation matrix has missing cell type pairs")
  if (max(abs(mat - t(mat))) > 1e-6) stop(panel, ": correlation matrix is not symmetric")
  if (min(mat) < -1 || max(mat) > 1) stop(panel, ": correlation outside [-1, 1]")
  message(panel, ": off-diagonal range = ", paste(signif(range(mat[upper.tri(mat)]), 3), collapse = " ~ "))

  ### Color dot annotation ####
  dot_colors <- celltype_colors[colnames(mat)]
  top_anno <- HeatmapAnnotation(
    show_annotation_name = FALSE,
    celltype = anno_points(x = rep(0.5, 7), gp = gpar(col = dot_colors, pch = 16), size = unit(10, "mm"),
                           ylim = c(0, 1), axis = FALSE, border = FALSE)
  )
  left_anno <- rowAnnotation(
    show_annotation_name = FALSE,
    celltype = anno_points(x = rep(0.5, 7), gp = gpar(col = dot_colors, pch = 16), size = unit(10, "mm"),
                           ylim = c(0, 1), axis = FALSE, border = FALSE)
  )
  celltype_legend <- Legend(
    labels = celltype_order, title = "Cell Type", type = "points", pch = 16, size = unit(5, "mm"),
    legend_gp = gpar(col = celltype_colors[celltype_order]),
    labels_gp = gpar(fontsize = 13), title_gp = gpar(fontsize = 15, fontface = "plain"),
    grid_height = unit(5, "mm"), grid_width = unit(5, "mm"), gap = unit(3, "mm"), background = "white"
  )
  ht <- Heatmap(
    mat,
    name = "Correlation",
    column_title = "All Cancer Types",
    top_annotation = top_anno,
    left_annotation = left_anno,
    col = col_fun,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title_gp = gpar(fontsize = 20, fontface = "bold"),
    heatmap_legend_param = list(
      direction = "horizontal", title_position = "topcenter", title_gp = gpar(fontface = "plain"),
      at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5)
    )
  )

  pdf(file.path(panels[[panel]]$plot_dir, paste0(panel, "_global_correlation.pdf")), width = 5, height = 4.5)
  draw(ht, show_heatmap_legend = TRUE, show_annotation_legend = TRUE,
       annotation_legend_list = list(celltype_legend), heatmap_legend_side = "bottom")
  decorate_heatmap_body("Correlation", {
    grid.rect(gp = gpar(col = "black", lwd = 2, fill = NA))
  })
  dev.off()
}
message("Fig.3a / Ext.Fig.7a complete: ", output_dir)
