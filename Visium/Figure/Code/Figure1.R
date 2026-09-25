# 분석 목적
#   - Source_Data_Fig1.xlsx(Fig.1a)로부터 암종별 scRNA-seq/Visium 샘플·세포·spot 수를 재현한다.
#
# 분석 흐름
#   1. Fig.1a 시트를 읽어 암종별 표본 규모 표를 정리한다.
#   2. scRNA-seq(샘플 수/세포 수)와 Visium(샘플 수/spot 수) 막대그래프를 각각 그린다.
#   3. 네 지표를 하나의 patchwork 그림으로 합쳐 저장한다.
#
# 주요 출력
#   - Fig1a.pdf (scRNA-seq/Visium 샘플·세포·spot 수 막대그래프)
#   - Fig1a_sample_summary.csv
#
# 출력 위치
#   - Output/Figure1/
#
# 참고: 원본 Fig.1a는 BioRender로 제작한 인체 실루엣을 중심에 둔 원형(fan) 막대그래프이며,
# Fig.1b는 전체가 BioRender 도식(workflow schematic)이다. 저작권이 있는 BioRender 아트워크와
# 수작업 배치는 코드로 재현하지 않고, 여기서는 Fig.1a에 담긴 수치 데이터만 동일한 정보량의
# 막대그래프로 재현한다. Fig.1b는 재현 대상 데이터가 없어 제외한다.

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
