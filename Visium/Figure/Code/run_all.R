# Paper_info/Source_Data의 공식 Source Data(Fig.3a는 번들 CSV)를 입력으로 사용해 Figure1, 2, 3(a), 4(c), 5(a,b)와
# Ext.Fig.7a의 R 패널을 순서대로 실행한다. Fig.3(b-h)와 Fig.4(a,b,d,f,g)/Fig.5(c,d)는 Python으로
# 재현되며 run_all_python.py가 담당한다 (겹치는 패널은 그쪽이 우선; README 참고).
# 개별 스크립트도 `Rscript Code/FigureN.R` 또는 `Rscript Code/FigureN/*.R`로 독립 실행 가능하다.

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
