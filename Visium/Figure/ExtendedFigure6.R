# 분석 목적
#   - Source_Data_Extended_Fig6.xlsx로부터 Extended Fig.6(top/bottom)을 재현한다
#     (Malignant/Normal compartment의 co-enrichment network).
#
# 분석 흐름
#   1. Malignant/Normal 각각의 edge list(이미 significance-filtered)를 불러온다.
#   2. 노드를 lineage로 색칠한 force-directed network를 그린다.
#
# 주요 출력
#   - ExtFig6_Malignant.pdf, ExtFig6_Normal.pdf
#
# 출력 위치
#   - Output/ExtendedFigure6/
#
# 참고: 원본은 Cytoscape에서 수작업으로 배치한 다이어그램이라 노드 좌표는 동일하지
# 않다. Source Data의 edge list로 force-directed layout(ggraph, layout="fr")을 새로
# 계산했으므로 노드 좌표는 원본과 다르지만 네트워크 구조(연결 관계)는 동일하다.

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(igraph); library(ggraph); library(ggplot2)
})
script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(script_path), "common/config.R"))

out_dir  <- file.path(output_dir, "ExtendedFigure6/Tables")
plot_dir <- file.path(output_dir, "ExtendedFigure6/Plots")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

subtype_lineage <- load_subtype_lineage()

plot_network <- function(sheet, panel, title) {
  # 시트 1행에 compartment 이름("Malignant"/"Normal")이 있고 2행이 실제 헤더다.
  dat <- read_source("Source_Data_Extended_Fig6.xlsx", sheet)
  names(dat) <- as.character(unlist(dat[1, ]))
  dat <- dat[-1, ]
  dat <- dat[!is.na(dat$node_a) & !is.na(dat$node_b), ]
  dat$pan_cancer_median_jaccard <- as.numeric(dat$pan_cancer_median_jaccard)
  write.csv(dat, file.path(out_dir, paste0(panel, "_network_edges.csv")), row.names = FALSE)

  g <- graph_from_data_frame(dat[, c("node_a", "node_b", "pan_cancer_median_jaccard")], directed = FALSE)
  V(g)$lineage <- subtype_lineage[V(g)$name]
  set.seed(1)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = pan_cancer_median_jaccard), color = "grey70", alpha = 0.4) +
    scale_edge_width(range = c(0.1, 1.2), guide = "none") +
    geom_node_point(aes(color = lineage), size = 3) +
    geom_node_text(aes(label = name), size = 2, repel = TRUE, max.overlaps = 15) +
    scale_color_manual(values = celltype_colors_spaced, name = "Lineage") +
    theme_void() +
    labs(title = title)
  # subtype 이름에 그리스 문자(gamma-delta T)가 있어 cairo_pdf 장치를 사용한다.
  cairo_pdf(file.path(plot_dir, paste0(panel, ".pdf")), width = 12, height = 12)
  print(p)
  dev.off()
}

plot_network("ED Fig.6 top", "ExtFig6_Malignant", "Malignant")
plot_network("ED Fig.6 bottom", "ExtFig6_Normal", "Normal")

message("Extended Figure 6 complete: ", dirname(out_dir))
