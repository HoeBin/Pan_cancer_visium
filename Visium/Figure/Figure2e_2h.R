# 분석 목적
#   - plot-ready proportion CSV에서 Figure 2e와 2h를 재현한다.
#
# 분석 흐름
#   1. scRNA-seq·Visium sample별 cell-type proportion을 불러온다.
#   2. 암종별 boxplot을 생성하고 요약 통계를 CSV로 저장한다.
#   3. compartment별 global/myeloid median proportion을 min–max scaling한다.
#   4. Figure 2h heatmap과 정규화 전·후 CSV를 저장한다.
#
# 주요 출력
#   - Fig2e_scRNAseq, Fig2e_Visium, Fig2h PDF
#   - 암종·cell type별 proportion summary 및 heatmap values CSV
#
# 출력 위치
#   - Figure/Output/Figure2/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2)
  library(ComplexHeatmap); library(circlize); library(viridisLite); library(grid); library(gridExtra)
})
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg)); root_dir <- dirname(script_path)
source(file.path(root_dir, "common/config.R"))
# "Load visium intermediate data" ---------------------------------------------
load_visium_intermediate <- function(include_subtypes = TRUE) {
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
  })

  obs_root <- DECON_DIR
  malignancy_file <- file.path(PROCESSED_DATA_DIR, "region_output_2025_02_02.csv")
  organ_map_file <- file.path(PROJECT_DIR, "Cell2Location_SubNum4/OrganToType.txt")

  obs_files <- list.files(
    obs_root,
    pattern = "Obs[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  obs_files <- obs_files[grepl("cell2location_map_.*_30000epoch/Obs[.]csv$", obs_files)]
  if (length(obs_files) != 15) {
    stop("Expected 15 intermediate Obs.csv files, found ", length(obs_files))
  }

  header <- names(data.table::fread(obs_files[1], nrows = 0))
  global_cols <- c(
    "Epi_enriched", "T_enriched", "B_enriched", "Mye_enriched",
    "EC_enriched", "Fib_enriched", "Mu_enriched"
  )
  first_subtype <- match("Age-associated_B", header)
  last_subtype <- match("vCAF", header)
  subtype_cols <- header[seq.int(first_subtype, last_subtype)]
  selected <- c(
    "V1", "cancer_type", "sample_id", "Data_GEO", "Sample_GEO",
    if (include_subtypes) subtype_cols else character(), global_cols
  )

  obs <- data.table::rbindlist(lapply(obs_files, function(path) {
    data.table::fread(path, select = selected, showProgress = FALSE)
  }), use.names = TRUE, fill = TRUE)
  data.table::setnames(obs, "V1", "cellid")

  organ_map <- data.table::fread(organ_map_file)
  obs[, Organ := cancer_type]
  obs <- merge(obs, organ_map, by = "Organ", all.x = TRUE)
  obs[!is.na(Cancer_type), cancer_type := Cancer_type]
  obs <- obs[cancer_type != "ESCA"]
  obs[, Cancer_type := NULL]

  malignancy <- data.table::fread(malignancy_file)
  id_col <- names(malignancy)[1]
  data.table::setnames(malignancy, id_col, "cellid")
  malignancy <- malignancy[, .(cellid, Malignancy = LocationNew, cnv_score)]
  obs <- merge(obs, malignancy, by = "cellid", all = FALSE)

  obs <- obs[
    sample_id != "GSM7757981_R3-S1" &
      !Sample_GEO %in% c(
        "SN123_A798015_Rep1", "SN124_A798015_Rep2",
        "GSE251950_21_01252_LI_SING"
      )
  ]
  obs[, Malignancy := factor(Malignancy, levels = c("Mal", "Bdy", "Normal"))]
  obs[]
}

global_celltype_columns <- function() {
  c(
    Epithelial = "Epi_enriched", Tcell = "T_enriched",
    Bcell = "B_enriched", Myeloid = "Mye_enriched",
    Endothelial = "EC_enriched", Fibroblast = "Fib_enriched",
    Mural = "Mu_enriched"
  )
}
out_dir <- file.path(root_dir, "Output/Figure2/Tables"); plot_dir <- file.path(root_dir, "Output/Figure2/Plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
base <- file.path(PROJECT_DIR, "Cell2Location_SubNum4/Supple_Table")
colors <- c(Epithelial = "red", Tcell = "blue", Bcell = "green", Myeloid = "purple", Endothelial = "magenta", Fibroblast = "orange", Mural = "cyan")
celltype_order <- c("Epithelial", "Tcell", "Bcell", "Myeloid", "Endothelial", "Fibroblast", "Mural")
cancer_order <- c("BRCA", "COCA", "HNCA", "KICA", "LICA", "LUCA", "OVCA", "PACA", "PRCA", "SKCA", "STCA", "THCA", "UECA")
celltype_labels <- c(Epithelial = "Epithelial", Tcell = "T cell", Bcell = "B cell", Myeloid = "Myeloid", Endothelial = "Endothelial", Fibroblast = "Fibroblast", Mural = "Mural")

# "Figure 2e" -----------------------------------------------------------------
sc <- fread(file.path(base, "Fig2e_sc.csv")); setnames(sc, c("Cancer_orig", "global_anno_edit"), c("cancer_type", "celltype"))
sc$celltype <- recode(sc$celltype, `T cell` = "Tcell", `B cell` = "Bcell")
vis <- fread(file.path(base, "Fig2e_visium.csv")); vis$celltype <- sub("_enriched$", "", vis$celltype)
vis$celltype <- recode(vis$celltype, Epi = "Epithelial", T = "Tcell", B = "Bcell", Mye = "Myeloid", EC = "Endothelial", Fib = "Fibroblast", Mu = "Mural")
sc <- sc[cancer_type != "ESCA"]; vis <- vis[cancer_type != "ESCA"]
sc$celltype <- factor(sc$celltype, levels = celltype_order)
vis$celltype <- factor(vis$celltype, levels = celltype_order)
sc$cancer_type <- factor(sc$cancer_type, levels = cancer_order)
vis$cancer_type <- factor(vis$cancer_type, levels = cancer_order)
summarise_prop <- function(x, platform) as.data.frame(x) %>% group_by(cancer_type, celltype) %>% summarise(platform = platform, n_samples = n(), median_proportion = median(proportion, na.rm = TRUE), mean_proportion = mean(proportion, na.rm = TRUE), q1 = quantile(proportion, 0.25, na.rm = TRUE), q3 = quantile(proportion, 0.75, na.rm = TRUE), .groups = "drop")
write.csv(bind_rows(summarise_prop(sc, "scRNA-seq"), summarise_prop(vis, "Visium")), file.path(out_dir, "Fig2e_key_results_proportion_by_cancer.csv"), row.names = FALSE)
write.csv(sc, file.path(out_dir, "Fig2e_scRNAseq_plot_data.csv"), row.names = FALSE); write.csv(vis, file.path(out_dir, "Fig2e_Visium_plot_data.csv"), row.names = FALSE)

plot_prop <- function(dat, title) ggplot(dat, aes(celltype, proportion)) + geom_jitter(size = 0.5) + geom_boxplot(aes(fill = celltype), outlier.shape = NA, alpha = 0.8) + facet_wrap(~cancer_type, ncol = 5, scales = "free_y", drop = TRUE) + scale_fill_manual(values = colors, drop = FALSE) + scale_x_discrete(labels = celltype_labels, drop = FALSE) + theme_bw() + labs(title = title, x = NULL, y = NULL) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 15, face = "plain", color = "black"), axis.title = element_text(size = 15, color = "black"), axis.text = element_text(size = 10, color = "black"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text = element_text(size = 10, face = "plain"))
pdf(file.path(plot_dir, "Fig2e_scRNAseq.pdf"), width = 7, height = 4); print(plot_prop(sc, "scRNA-seq")); dev.off()
pdf(file.path(plot_dir, "Fig2e_Visium.pdf"), width = 7, height = 4); print(plot_prop(vis, "Visium")); dev.off()

# "Figure 2h" -----------------------------------------------------------------
# 원본 생성 스크립트(Malignancy_Analysis/Proportion_heatmap.R)와 동일하게, 캐시된 요약
# CSV가 아니라 원시 spot 단위 중간 데이터를 Malignancy별로 직접 pooling해 median을 낸다.
obs_h <- load_visium_intermediate(include_subtypes = TRUE)
myeloid_subtypes <- c("CXCL3+_Macro", "CD14+_Mono", "CD1C+_cDC", "FOLR2+_Macro", "HSPA6+_Macro",
                       "SPP1+_cycMacro", "CD16+_Mono", "Macro", "LILRA4+_pDC", "CLEC9A+_cDC",
                       "TNFSF10+_Macro", "SPP1+_CXCL3+_Macro", "SPP1+_Macro", "SPP1+_MT1+_Macro",
                       "LAMP3+_cDC", "Mast")
type_columns <- list(Global = global_celltype_columns(), Myeloid = setNames(myeloid_subtypes, myeloid_subtypes))
scale01 <- function(x) if (diff(range(x, na.rm = TRUE)) == 0) rep(0, length(x)) else (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))
mal_compartment <- c(Mal = "Malignant", Bdy = "Boundary", Normal = "Normal")

h <- bind_rows(lapply(names(type_columns), function(type) {
  cols <- type_columns[[type]]
  summary_wide <- as.data.frame(obs_h) %>% group_by(Malignancy) %>%
    summarise(across(all_of(unname(cols)), ~ median(.x, na.rm = TRUE)), .groups = "drop")
  long <- pivot_longer(summary_wide, -Malignancy, names_to = "column", values_to = "Median")
  long$CellType <- names(cols)[match(long$column, cols)]
  long$Type <- type
  long$Compartment <- mal_compartment[as.character(long$Malignancy)]
  long %>% group_by(CellType) %>% mutate(ScaledMedian = scale01(Median)) %>% ungroup() %>%
    select(Type, Malignancy, Compartment, CellType, Median, ScaledMedian)
}))
write.csv(h, file.path(out_dir, "Fig2h_key_results_compartment_proportions.csv"), row.names = FALSE)

col_fun <- colorRamp2(c(0, 0.9, 1), viridis(3, option = "C"))
heatmaps <- list()
for (type in names(type_columns)) {
  x <- filter(h, Type == type) %>% select(CellType, Compartment, ScaledMedian) %>% pivot_wider(names_from = Compartment, values_from = ScaledMedian) %>% as.data.frame()
  rownames(x) <- x$CellType; x$CellType <- NULL
  mat <- as.matrix(x[, c("Malignant", "Boundary", "Normal"), drop = FALSE])
  heatmaps[[type]] <- Heatmap(mat, name = "proportion", column_title = type,
    col = col_fun, cluster_columns = FALSE, cluster_rows = TRUE,
    column_names_rot = 45, row_names_gp = gpar(fontsize = 15),
    column_names_gp = gpar(fontsize = 15),
    column_title_gp = gpar(fontsize = 20, fontface = "bold"),
    heatmap_legend_param = list(at = c(0, 1), labels = c("min", "max"),
      title_gp = gpar(fontsize = 15, fontface = "bold"),
      labels_gp = gpar(fontsize = 12)))
}
pdf(file.path(plot_dir, "Fig2h.pdf"), width = 7, height = 7)
for (type in names(heatmaps)) {
  draw(heatmaps[[type]], padding = unit(c(2, 40, 2, 2), "mm"))
}
dev.off()
message("Figure 2e and 2h complete: ", dirname(out_dir))
