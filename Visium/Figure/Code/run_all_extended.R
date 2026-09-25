# Runs the reproducible R panels of Extended Data Figures 1, 2, 4 and 5(a) in order, using only the officially
# distributed Source Data in Paper_info/Source_Data as input.
# Ext. Fig. 3 (companion of Fig. 2), Ext. Fig. 5b and Ext. Fig. 6-10 are reproduced in Python and handled by
# run_all_python.py (Python takes precedence for overlapping panels; see the README).
# Each script can also be run on its own with `Rscript Code/ExtendedFigureN.R`.

script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
root_dir    <- dirname(script_path)

figure_scripts <- file.path(root_dir, paste0("ExtendedFigure", c(1, 2, 4, 5), ".R"))
for (script in figure_scripts) {
  message("==== Running ", basename(script), " ====")
  system2("Rscript", shQuote(script))
}
