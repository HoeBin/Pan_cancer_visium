# Paper_info/Source_Data의 공식 Source Data만을 입력으로 사용해 Extended Data
# Figure 1, 2, 4, 5(a)의 재현 가능한 R 패널을 순서대로 실행한다.
# Ext.Fig.3(companion of Fig.2), Ext.Fig.5b, Ext.Fig.6-10은 Python으로 재현되며
# run_all_python.py가 담당한다 (겹치는 패널은 그쪽이 우선; README 참고).
# 개별 스크립트도 `Rscript Code/ExtendedFigureN.R`로 독립 실행 가능하다.

script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
root_dir    <- dirname(script_path)

figure_scripts <- file.path(root_dir, paste0("ExtendedFigure", c(1, 2, 4, 5), ".R"))
for (script in figure_scripts) {
  message("==== Running ", basename(script), " ====")
  system2("Rscript", shQuote(script))
}
