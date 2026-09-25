# Purpose
#   - Reproduces Extended Fig.2a-g from Source_Data_Extended_Fig2.xlsx
#     (subtype x cancer type Ro/e (observed/expected ratio) heatmap per major cell type).
#
# Workflow
#   1. Build a Cell_subtype x Cancer_type Ro/e matrix for each of the 7 Major_cell_types.
#   2. Draw a heatmap that shows the Ro/e value and significance symbol (Ro_e_symbol) in each cell.
#
# Main outputs
#   - ExtFig2a(T cell) … ExtFig2g(Epithelial) PDFs (7 files)
#
# Output location
#   - Output/ExtendedFigure2/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ComplexHeatmap); library(circlize)
})
script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(script_path), "common/config.R"))

out_dir  <- file.path(output_dir, "ExtendedFigure2/Tables")
plot_dir <- file.path(output_dir, "ExtendedFigure2/Plots")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

ed2 <- read_source("Source_Data_Extended_Fig2.xlsx", "ED Fig.2")
ed2 <- ed2[!is.na(ed2$Cell_subtype) & ed2$Cancer_type %in% cancer_order, ]
write.csv(ed2, file.path(out_dir, "ExtFig2_Roe_by_cancer.csv"), row.names = FALSE)

panels <- c("a" = "T cell", "b" = "B cell", "c" = "Endothelial", "d" = "Myeloid",
            "e" = "Fibroblast", "f" = "Mural", "g" = "Epithelial")

for (panel in names(panels)) {
  major <- panels[panel]
  dat <- ed2[ed2$Major_cell_type == major, ]
  subtype_order <- sort(unique(dat$Cell_subtype))
  dat$Cell_subtype <- factor(dat$Cell_subtype, levels = subtype_order)
  dat$Cancer_type  <- factor(dat$Cancer_type, levels = cancer_order)

  mat   <- xtabs(Ro_e ~ Cell_subtype + Cancer_type, data = dat)
  label_dat <- dat %>% mutate(cell_label = paste0(Ro_e_display, "\n", Ro_e_symbol))
  label_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  for (i in seq_len(nrow(label_dat))) {
    label_mat[as.character(label_dat$Cell_subtype[i]), as.character(label_dat$Cancer_type[i])] <- label_dat$cell_label[i]
  }

  col_fun <- colorRamp2(c(0, 1, max(mat, na.rm = TRUE)), c("blue", "white", "red"))
  ht <- Heatmap(mat, name = "Ro/e", col = col_fun,
    column_title = major,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), column_names_rot = 45,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(label_mat[i, j], x, y, gp = gpar(fontsize = 6.5, col = "black"))
    })
  # Use cairo_pdf: subtype names contain Greek letters (e.g. gamma-delta T), which the default pdf()
  # device (mbcsToSbcs) renders incorrectly.
  cairo_pdf(file.path(plot_dir, paste0("ExtFig2", panel, ".pdf")),
      width = 9, height = max(4, 0.35 * length(subtype_order) + 1.5))
  draw(ht)
  dev.off()
}

message("Extended Figure 2 complete: ", dirname(out_dir))
