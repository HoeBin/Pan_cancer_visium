# Figure 스크립트 공용 설정: 경로, 색상 팔레트, Source Data 로더.
# 각 FigureN.R은 source(file.path(root_dir, "common/config.R"))로 이 파일을 불러온다.
# 논문에 공식 배포된 Source Data(Figure/Source_Data/*.xlsx)만 입력으로 사용하는
# self-contained 재현 스크립트다. 저장소를 clone한 뒤 별도 경로 설정 없이 바로
# 실행할 수 있도록 모든 경로를 이 저장소 안(Figure/) 기준 상대경로로 계산한다.

script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
root_dir        <- dirname(script_path)                    # .../Visium/Figure
source_data_dir <- file.path(root_dir, "Source_Data")
output_dir      <- file.path(root_dir, "Output")

# "Source Data 로더" -----------------------------------------------------------
read_source <- function(fig_file, sheet) {
  suppressPackageStartupMessages(library(readxl))
  path <- file.path(source_data_dir, fig_file)
  as.data.frame(read_excel(path, sheet = sheet))
}

# "색상 팔레트 및 순서" ----------------------------------------------------------
celltype_order  <- c("Epithelial", "Tcell", "Bcell", "Myeloid", "Endothelial", "Fibroblast", "Mural")
celltype_labels <- c(Epithelial = "Epithelial", Tcell = "T cell", Bcell = "B cell", Myeloid = "Myeloid",
                      Endothelial = "Endothelial", Fibroblast = "Fibroblast", Mural = "Mural")
celltype_colors <- c(Epithelial = "#FF0000", Tcell = "#0000FF", Bcell = "#00FF00", Myeloid = "#9932CC",
                      Endothelial = "#FF00FF", Fibroblast = "#FFA500", Mural = "#00EEEE")
# Source Data의 "T cell"/"B cell" 표기(공백 포함)에 대응하는 별칭
celltype_colors_spaced <- setNames(celltype_colors, celltype_labels[names(celltype_colors)])

mal_order  <- c("Malignant", "Boundary", "Normal")
mal_colors <- c(Malignant = "#FF0000", Boundary = "#008B00", Normal = "#4169E1")

cancer_order <- c("BRCA", "COCA", "HNCA", "KICA", "LICA", "LUCA", "OVCA", "PACA", "PRCA", "SKCA", "STCA", "THCA", "UECA")

mal_comparisons <- list(c("Normal", "Boundary"), c("Boundary", "Malignant"), c("Normal", "Malignant"))

# "Cell subtype -> lineage 매핑" -----------------------------------------------
# ED Fig.2(Major_cell_type/Cell_subtype)에서 유도한 93개 세부 subtype -> 7개 lineage
# 매핑. Extended Fig4d(상관관계 bar plot)와 Fig6(co-enrichment network)의 노드 색상에
# 공통으로 사용한다. subtype 이름 표기가 스크립트마다 "_"/" " 혼용이라 둘 다 키로 등록한다.
load_subtype_lineage <- function() {
  ed2 <- read_source("Source_Data_Extended_Fig2.xlsx", "ED Fig.2")
  ed2 <- ed2[!is.na(ed2$Cell_subtype), c("Major_cell_type", "Cell_subtype")]
  ed2 <- unique(ed2)
  lookup <- setNames(ed2$Major_cell_type, ed2$Cell_subtype)
  lookup_us <- setNames(ed2$Major_cell_type, gsub(" ", "_", ed2$Cell_subtype, fixed = TRUE))
  lookup <- c(lookup, lookup_us[!names(lookup_us) %in% names(lookup)])
  # gamma-delta T의 ASCII 표기("gd_T")를 쓰는 시트가 있어 별칭으로 추가한다.
  c(lookup, "gd_T" = "T cell")
}

# "공통 ggplot theme" ------------------------------------------------------------
base_theme <- ggplot2::theme_bw() + ggplot2::theme(
  strip.background = ggplot2::element_rect(fill = "lightgray", color = "black", linewidth = 0.5),
  strip.text       = ggplot2::element_text(size = 15, face = "plain", color = "black"),
  axis.text        = ggplot2::element_text(size = 15, color = "black"),
  axis.title       = ggplot2::element_text(size = 15),
  plot.title       = ggplot2::element_text(hjust = 0.5, size = 15, face = "plain"),
  legend.position  = "right"
)
