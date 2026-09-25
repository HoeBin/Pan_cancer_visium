# Runs the R panels (Figures 1, 2, 3a, 4c, 5a-b and Ext. Fig. 7a) in order, using the officially distributed
# Source Data in Paper_info/Source_Data (Fig. 3a uses a bundled CSV) as input. Fig. 3(b-h), Fig. 4(a,b,d,f,g) and
# Fig. 5(c,d) are reproduced in Python and handled by run_all_python.py (Python takes precedence for overlapping
# panels; see the README).
# Each script can also be run on its own with `Rscript Code/FigureN.R` or `Rscript Code/FigureN/*.R`.

script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
root_dir    <- dirname(script_path)

figure_scripts <- c(
  file.path(root_dir, "Figure1.R"),
  file.path(root_dir, "Figure2.R"),
  file.path(root_dir, "Figure3/Fig3a_ExtDataFig7a_correlation_heatmap.R"),
  file.path(root_dir, "Figure4/Fig4c_pEMT_abundance_correlation.R"),
  file.path(root_dir, "Figure5/Fig5ab_pancancer_LR_zscore.R")
)
for (script in figure_scripts) {
  message("==== Running ", basename(script), " ====")
  system2("Rscript", shQuote(script))
}
