# Shared settings for the figure scripts: paths, color palettes and the Source Data loader.
# FigureN.R directly under Code/ loads this file with source(file.path(root_dir, "common/config.R")).
# Scripts one level deeper, such as Code/FigureN/0X_*.R, define code_dir (= Code/) themselves before sourcing —
# commandArgs() only reports the path of the top-level script that was run, so at a different depth the
# dirname() calculation below would not match the script location.

script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
if (!exists("code_dir")) code_dir <- dirname(script_path)  # .../Github_code_v2/Code
root_dir       <- code_dir
github_code_dir <- dirname(code_dir)                        # .../Github_code_v2
source_data_dir <- file.path(github_code_dir, "Paper_info/Paper_info/Source_Data")
output_dir      <- file.path(github_code_dir, "Output")

# "Source Data loader" -----------------------------------------------------------
read_source <- function(fig_file, sheet) {
  suppressPackageStartupMessages(library(readxl))
  path <- file.path(source_data_dir, fig_file)
  as.data.frame(read_excel(path, sheet = sheet))
}

# "Color palettes and order" ----------------------------------------------------------
celltype_order  <- c("Epithelial", "Tcell", "Bcell", "Myeloid", "Endothelial", "Fibroblast", "Mural")
celltype_labels <- c(Epithelial = "Epithelial", Tcell = "T cell", Bcell = "B cell", Myeloid = "Myeloid",
                      Endothelial = "Endothelial", Fibroblast = "Fibroblast", Mural = "Mural")
celltype_colors <- c(Epithelial = "#FF0000", Tcell = "#0000FF", Bcell = "#00FF00", Myeloid = "#9932CC",
                      Endothelial = "#FF00FF", Fibroblast = "#FFA500", Mural = "#00EEEE")
# Alias for the "T cell"/"B cell" labels (with a space) used in Source Data
celltype_colors_spaced <- setNames(celltype_colors, celltype_labels[names(celltype_colors)])

mal_order  <- c("Malignant", "Boundary", "Normal")
mal_colors <- c(Malignant = "#FF0000", Boundary = "#008B00", Normal = "#4169E1")

cancer_order <- c("BRCA", "COCA", "HNCA", "KICA", "LICA", "LUCA", "OVCA", "PACA", "PRCA", "SKCA", "STCA", "THCA", "UECA")

mal_comparisons <- list(c("Normal", "Boundary"), c("Boundary", "Malignant"), c("Normal", "Malignant"))

# "Cell subtype -> lineage mapping" -----------------------------------------------
# Mapping of the 93 fine-grained subtypes to the 7 lineages, derived from ED Fig.2 (Major_cell_type/Cell_subtype).
# Shared by the node colors of Extended Fig4d (correlation bar plot) and Fig6 (co-enrichment network).
# Subtype names mix "_" and " " across scripts, so both spellings are registered as keys.
load_subtype_lineage <- function() {
  ed2 <- read_source("Source_Data_Extended_Fig2.xlsx", "ED Fig.2")
  ed2 <- ed2[!is.na(ed2$Cell_subtype), c("Major_cell_type", "Cell_subtype")]
  ed2 <- unique(ed2)
  lookup <- setNames(ed2$Major_cell_type, ed2$Cell_subtype)
  lookup_us <- setNames(ed2$Major_cell_type, gsub(" ", "_", ed2$Cell_subtype, fixed = TRUE))
  lookup <- c(lookup, lookup_us[!names(lookup_us) %in% names(lookup)])
  # Some sheets use the ASCII spelling of gamma-delta T ("gd_T"), so add it as an alias.
  c(lookup, "gd_T" = "T cell")
}

# "Common ggplot theme" ------------------------------------------------------------
base_theme <- ggplot2::theme_bw() + ggplot2::theme(
  strip.background = ggplot2::element_rect(fill = "lightgray", color = "black", linewidth = 0.5),
  strip.text       = ggplot2::element_text(size = 15, face = "plain", color = "black"),
  axis.text        = ggplot2::element_text(size = 15, color = "black"),
  axis.title       = ggplot2::element_text(size = 15),
  plot.title       = ggplot2::element_text(hjust = 0.5, size = 15, face = "plain"),
  legend.position  = "right"
)
