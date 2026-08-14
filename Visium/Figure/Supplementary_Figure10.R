# 분석 목적
#   - compartment별 cell subtype proportion으로 Supplementary Fig.10을 재현한다.
#
# 분석 흐름
#   1. Obs.csv와 malignancy 중간 결과를 결합한다.
#   2. subtype 및 암종별 median proportion을 계산한다.
#   3. 각 subtype/암종 내 min–max scaling을 적용한다.
#   4. heatmap PDF와 원값·정규화값 CSV를 저장한다.
#
# 주요 출력
#   - Supplementary_Fig10a,b PDF
#   - median proportion 및 min–max scaled proportion CSV
#
# 출력 위치
#   - Figure/Output/Supplementary_Figure10/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(viridisLite)
  library(grid)
  library(gridExtra)
})
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
root_dir <- dirname(script_path)
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
out_dir <- file.path(root_dir, "Output/Supplementary_Figure10/Tables")
plot_dir <- file.path(root_dir, "Output/Supplementary_Figure10/Plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# "Load and harmonize intermediate data" -------------------------------------
obs <- load_visium_intermediate(include_subtypes = TRUE)
celltype_map <- fread(file.path(PROJECT_DIR, "Cell2Location_SubNum4/Data_Info/celltype_mapping.csv"))
rename_map <- setNames(celltype_map$after_celltype, celltype_map$before_celltype)
matched <- intersect(names(rename_map), names(obs))
setnames(obs, matched, unname(rename_map[matched]))

groups <- list(
  Epithelial = paste0("Epi", 0:16),
  Tcell = c("Trm_cytotoxic_T", "Tem/Trm_cytotoxic_T", "CD16-_NK", "exT", "Tem/Effector_helper_T", "regT", "Tem/Temra_cytotoxic_T", "ILC3", "Follicular_helper_T", "Tcm/Naive_helper_T", "FOXP3-_regT", "MAIT", "gamma-delta_T", "Type_17_helper_T", "NK", "Tcm/Naive_cytotoxic_T", "Regulatory_T", "cycling_exT", "CD16+_NK", "ISG15+_T", "cycling_regT", "Type_1_helper_T"),
  Bcell = c("B", "Naive_B", "Memory_B", "Germinal_center_B", "Plasma", "Proliferative_germinal_center_B", "Age-associated_B", "Plasmablasts"),
  Myeloid = c("CXCL3+_Macro", "CD14+_Mono", "CD1C+_cDC", "FOLR2+_Macro", "HSPA6+_Macro", "SPP1+_cycMacro", "CD16+_Mono", "Macro", "LILRA4+_pDC", "CLEC9A+_cDC", "TNFSF10+_Macro", "SPP1+_CXCL3+_Macro", "SPP1+_Macro", "SPP1+_MT1+_Macro", "LAMP3+_cDC", "Mast"),
  Endothelial = c("Venous_EC", "Capillary_EC", "Arterial_EC", "Fibrosis_PGF+_Tip_EC", "Venous_iEC", "PGF+_Tip_EC", "Lymphatic_EC", "HMGB2+_cycEC", "Tip_EC", "SPRY1+_EC", "CCL2+_iEC", "HMOX1+_iEC", "CD14+_circEC", "TMEM100+_EC"),
  Fibroblast = c("CXCL14+_mCAF", "PI16+_iCAF", "iCAF", "Fibroblast", "apCAF", "vCAF", "IL6+_iCAF", "pnCAF", "mCAF", "ISG15+_mCAF", "tCAF", "HSP+_tCAF", "cyc_mCAF", "myoCAF"),
  Mural = c("SMC", "Pericyte")
)
groups <- lapply(groups, function(x) {
  mapped <- ifelse(x %in% names(rename_map), unname(rename_map[x]), x)
  intersect(mapped, names(obs))
})
compartment_names <- c(Mal = "Malignant", Bdy = "Boundary", Normal = "Normal")
scale01 <- function(x) if (all(is.na(x)) || diff(range(x, na.rm = TRUE)) == 0) rep(0, length(x)) else (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))
col_fun <- colorRamp2(c(0, 0.5, 1), viridis(3, option = "C"))
grDevices::pdf(NULL)

capture_heatmap <- function(mat, title, row_font = 8) {
  grid.grabExpr(draw(Heatmap(mat, column_title = title, name = "proportion", col = col_fun,
    show_row_dend = TRUE, cluster_columns = FALSE, show_column_dend = TRUE, cluster_rows = TRUE,
    row_names_gp = gpar(fontsize = row_font), column_names_gp = gpar(fontsize = 10),
    column_title_gp = gpar(fontsize = 13), column_names_rot = 45,
    heatmap_legend_param = list(title_gp = gpar(fontsize = 9))),
    padding = unit(c(2, 20, 2, 2), "mm")),
    device = function(width, height) grDevices::pdf(NULL, width = width, height = height))
}

# "Panel 10a: subtype heatmaps" -----------------------------------------------
panel10a_groups <- c("Epithelial", "Tcell", "Bcell", "Endothelial", "Fibroblast", "Mural")
raw_a <- list(); scaled_a <- list(); grobs_a <- list()
for (group_name in panel10a_groups) {
  columns <- groups[[group_name]]
  summary_long <- as.data.frame(obs)[, c("Malignancy", columns)] %>%
    group_by(Malignancy) %>% summarise(across(all_of(columns), median, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(-Malignancy, names_to = "subtype", values_to = "median_proportion")
  summary_long$compartment <- unname(compartment_names[as.character(summary_long$Malignancy)])
  summary_long$global_celltype <- group_name
  summary_long <- summary_long %>% group_by(subtype) %>% mutate(scaled_proportion = scale01(median_proportion)) %>% ungroup()
  raw_a[[group_name]] <- select(summary_long, global_celltype, subtype, compartment, median_proportion)
  scaled_a[[group_name]] <- select(summary_long, global_celltype, subtype, compartment, scaled_proportion)
  mat <- summary_long %>% select(subtype, compartment, scaled_proportion) %>% pivot_wider(names_from = compartment, values_from = scaled_proportion) %>% as.data.frame()
  rownames(mat) <- mat$subtype; mat$subtype <- NULL
  mat <- as.matrix(mat[, c("Malignant", "Boundary", "Normal"), drop = FALSE])
  grobs_a[[group_name]] <- capture_heatmap(mat, group_name, ifelse(nrow(mat) > 15, 6, 8))
}
write.csv(bind_rows(raw_a), file.path(out_dir, "Supplementary_Fig10a_median_proportions.csv"), row.names = FALSE)
write.csv(bind_rows(scaled_a), file.path(out_dir, "Supplementary_Fig10a_scaled_proportions.csv"), row.names = FALSE)
pdf(file.path(plot_dir, "Supplementary_Fig10a.pdf"), width = 8, height = 7)
for (group_name in panel10a_groups) {
  grid.newpage()
  grid.draw(grobs_a[[group_name]])
}
dev.off()

# "Panel 10b: cancer-by-global-celltype heatmaps" -----------------------------
global_map <- global_celltype_columns()
global_display <- c(Epithelial = "Epithelial", Tcell = "Tcell", Bcell = "Bcell", Myeloid = "Myeloid", Endothelial = "Endothelial", Fibroblast = "Fibroblast", Mural = "Mural")
raw_b <- list(); scaled_b <- list(); grobs_b <- list()
for (group_name in names(global_map)) {
  column <- unname(global_map[group_name])
  summary_long <- as.data.frame(obs)[, c("cancer_type", "Malignancy", column)]
  names(summary_long)[3] <- "proportion"
  summary_long <- summary_long %>% group_by(cancer_type, Malignancy) %>% summarise(median_proportion = median(proportion, na.rm = TRUE), .groups = "drop") %>%
    filter(cancer_type != "ESCA") %>% mutate(compartment = unname(compartment_names[as.character(Malignancy)])) %>%
    group_by(cancer_type) %>% mutate(scaled_proportion = scale01(median_proportion)) %>% ungroup()
  summary_long$global_celltype <- group_name
  raw_b[[group_name]] <- select(summary_long, global_celltype, cancer_type, compartment, median_proportion)
  scaled_b[[group_name]] <- select(summary_long, global_celltype, cancer_type, compartment, scaled_proportion)
  mat <- summary_long %>% select(cancer_type, compartment, scaled_proportion) %>% pivot_wider(names_from = compartment, values_from = scaled_proportion) %>% as.data.frame()
  rownames(mat) <- mat$cancer_type; mat$cancer_type <- NULL
  mat <- as.matrix(mat[, c("Malignant", "Boundary", "Normal"), drop = FALSE])
  grobs_b[[group_name]] <- capture_heatmap(mat, global_display[group_name], 7)
}
write.csv(bind_rows(raw_b), file.path(out_dir, "Supplementary_Fig10b_median_proportions_nonESCA.csv"), row.names = FALSE)
write.csv(bind_rows(scaled_b), file.path(out_dir, "Supplementary_Fig10b_scaled_proportions_nonESCA.csv"), row.names = FALSE)
pdf(file.path(plot_dir, "Supplementary_Fig10b.pdf"), width = 7, height = 7)
for (group_name in names(global_map)) {
  grid.newpage()
  grid.draw(grobs_b[[group_name]])
}
dev.off()

message("Supplementary Figure 10 complete: ", dirname(out_dir))
