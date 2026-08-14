# 분석 목적
#   - Cell2location abundance 중간 결과에서 Supplementary Figures 14와 16을 재현한다.
#
# 분석 흐름
#   1. sample별 좌표와 global cell-type abundance를 불러온다.
#   2. k-NN(k=8) 공간 가중치로 Moran's I를 계산한다.
#   3. 전체 sample 및 compartment별 결과를 시각화한다.
#   4. Moran's I와 Wilcoxon 검정 결과를 CSV로 저장한다.
#
# 주요 출력
#   - Supplementary Fig.14a,b 및 Fig.16a,b PDF
#   - sample/compartment별 Moran's I와 통계 검정 CSV
#
# 출력 위치
#   - Figure/Output/Supplementary_Figure14_16/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(spdep)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
root_dir <- dirname(script_path)
source(file.path(root_dir, "common/config.R"))
out_dir  <- file.path(root_dir, "Output/Supplementary_Figure14_16/Tables")
plot_dir <- file.path(root_dir, "Output/Supplementary_Figure14_16/Plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

obs_root <- DECON_DIR
malignancy_file <- file.path(PROCESSED_DATA_DIR, "region_output_2025_02_02.csv")
organ_map_file <- file.path(PROJECT_DIR, "Cell2Location_SubNum4/OrganToType.txt")

# "Load abundance intermediates" ---------------------------------------------
obs_files <- list.files(obs_root, "Obs_abundance[.]csv$", recursive = TRUE, full.names = TRUE)
obs_files <- obs_files[grepl("cell2location_map_.*_30000epoch/Obs_abundance[.]csv$", obs_files)]
required <- c(
  "V1", "array_row", "array_col", "cancer_type", "sample_id",
  "Sample_GEO", "Epi_enriched", "T_enriched", "B_enriched",
  "Mye_enriched", "EC_enriched", "Fib_enriched", "Mu_enriched"
)
obs <- rbindlist(lapply(obs_files, function(path) {
  fread(path, select = required, showProgress = FALSE)
}), use.names = TRUE)
setnames(obs, c("V1", "Epi_enriched", "T_enriched", "B_enriched", "Mye_enriched",
                "EC_enriched", "Fib_enriched", "Mu_enriched"),
         c("cellid", "Epithelial_enriched", "Tcell_enriched", "Bcell_enriched",
           "Myeloid_enriched", "Endothelial_enriched", "Fibroblast_enriched", "Mural_enriched"))

organ_map <- fread(organ_map_file)
obs[, Organ := cancer_type]
obs <- merge(obs, organ_map, by = "Organ", all.x = TRUE)
obs[!is.na(Cancer_type), cancer_type := Cancer_type]
obs[, Cancer_type := NULL]
obs <- obs[cancer_type != "ESCA"]
obs <- obs[!Sample_GEO %in% c("SN123_A798015_Rep1", "SN124_A798015_Rep2", "GSE251950_21_01252_LI_SING")]

cell_cols <- paste0(c("Epithelial", "Tcell", "Bcell", "Myeloid", "Endothelial", "Fibroblast", "Mural"), "_enriched")
cell_levels <- c("Epithelial", "T cell", "B cell", "Myeloid", "Endothelial", "Fibroblast", "Mural")
cell_colors <- c(Epithelial = "red", `T cell` = "blue", `B cell` = "green", Myeloid = "purple",
                 Endothelial = "magenta", Fibroblast = "orange", Mural = "cyan")
organ_levels <- c("BRCA", "COCA", "KICA", "OVCA", "HNCA", "LUCA", "PACA", "STCA", "LICA", "PRCA", "SKCA", "THCA", "UECA")
organ_colors <- setNames(scales::hue_pal()(length(organ_levels)), organ_levels)

compute_morans <- function(df, columns, k = 8) {
  df <- df[complete.cases(df[, c("array_row", "array_col")]), , drop = FALSE]
  n <- nrow(df)
  if (n < 3) return(tibble(cell_type = columns, morans_I = NA_real_, p_value = NA_real_, n_points = n, k_used = NA_integer_))
  k_eff <- min(k, n - 1L)
  coords <- as.matrix(df[, c("array_row", "array_col")])
  lw <- nb2listw(knn2nb(knearneigh(coords, k = k_eff)), style = "W", zero.policy = TRUE)
  map_dfr(columns, function(column) {
    values <- df[[column]]
    if (anyNA(values) || length(unique(values)) < 2) return(tibble(cell_type = column, morans_I = NA_real_, p_value = NA_real_, n_points = n, k_used = k_eff))
    result <- moran.test(values, lw, randomisation = TRUE, alternative = "two.sided", zero.policy = TRUE)
    tibble(cell_type = column, morans_I = unname(result$estimate[["Moran I statistic"]]), p_value = result$p.value, n_points = n, k_used = k_eff)
  })
}

format_celltype <- function(x) {
  x <- sub("_enriched$", "", x)
  x[x == "Tcell"] <- "T cell"
  x[x == "Bcell"] <- "B cell"
  factor(x, levels = cell_levels)
}

# "Supplementary Figure 14" ---------------------------------------------------
filter_cache <- "--filter-cache" %in% commandArgs(trailingOnly = TRUE)
morans_sample_path <- file.path(out_dir, "Supplementary_Fig14_morans_by_sample.csv")
if (filter_cache && file.exists(morans_sample_path)) {
  morans_sample <- read.csv(morans_sample_path) %>% filter(cancer_type != "ESCA")
} else {
  morans_sample <- as.data.frame(obs) %>%
    group_by(cancer_type, sample_id) %>%
    group_modify(~ compute_morans(as.data.frame(.x), cell_cols, 8)) %>%
    ungroup()
}
morans_sample$cell_type <- format_celltype(morans_sample$cell_type)
write.csv(morans_sample, morans_sample_path, row.names = FALSE)

p14a <- ggplot(morans_sample, aes(cell_type, morans_I, fill = cell_type)) +
  geom_boxplot() + scale_fill_manual(values = cell_colors) + theme_cowplot() +
  labs(x = "cell_type", y = "morans_I") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 16), axis.title = element_text(size = 18))
p14b <- ggplot(filter(morans_sample, cancer_type != "ESCA"), aes(cancer_type, morans_I, color = cancer_type)) +
  geom_boxplot() + facet_wrap(~cell_type) + scale_color_manual(values = organ_colors) +
  theme_cowplot() + labs(x = NULL, y = "morans_I") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 7), strip.text = element_text(size = 9))

pdf(file.path(plot_dir, "Supplementary_Fig14a.pdf"), width = 4, height = 3.5); print(p14a); dev.off()
pdf(file.path(plot_dir, "Supplementary_Fig14b.pdf"), width = 6, height = 3); print(p14b); dev.off()

# "Supplementary Figure 16" ---------------------------------------------------
malignancy <- fread(malignancy_file)
setnames(malignancy, names(malignancy)[1], "cellid")
malignancy <- malignancy[, .(cellid, Malignancy = LocationNew)]
obs_mal <- merge(obs, malignancy, by = "cellid", all = FALSE)
obs_mal <- obs_mal[sample_id != "GSM7757981_R3-S1"]

morans_compartment_path <- file.path(out_dir, "Supplementary_Fig16_morans_by_compartment.csv")
if (filter_cache && file.exists(morans_compartment_path)) {
  morans_compartment <- read.csv(morans_compartment_path) %>% filter(cancer_type != "ESCA")
} else {
  morans_compartment <- as.data.frame(obs_mal) %>%
    group_by(cancer_type, sample_id, Malignancy) %>%
    group_modify(~ compute_morans(as.data.frame(.x), cell_cols, 8)) %>%
    ungroup() %>%
    filter(!is.na(p_value), p_value < 0.05)
}
morans_compartment$cell_type <- format_celltype(morans_compartment$cell_type)
morans_compartment$Malignancy <- factor(morans_compartment$Malignancy, levels = c("Mal", "Bdy", "Normal"))
write.csv(morans_compartment, morans_compartment_path, row.names = FALSE)

comparisons <- list(c("Normal", "Bdy"), c("Bdy", "Mal"), c("Normal", "Mal"))
test_all <- compare_means(morans_I ~ Malignancy, data = morans_compartment, method = "wilcox.test", comparisons = comparisons)
test_celltype <- compare_means(morans_I ~ Malignancy, data = morans_compartment, group.by = "cell_type", method = "wilcox.test", comparisons = comparisons)
write.csv(test_all, file.path(out_dir, "Supplementary_Fig16a_wilcoxon.csv"), row.names = FALSE)
write.csv(test_celltype, file.path(out_dir, "Supplementary_Fig16b_wilcoxon_by_celltype.csv"), row.names = FALSE)

mal_colors <- c(Mal = "#FF0000", Bdy = "#008B00", Normal = "#4169E1")
mal_labels <- c(Mal = "Malignant", Bdy = "Boundary", Normal = "Normal")
p16a <- ggplot(morans_compartment, aes(Malignancy, morans_I, color = Malignancy)) +
  geom_boxplot() + stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", label.y = 0.7, step.increase = c(0, 0.05, 0.07), tip.length = 0.01, size = 6) +
  scale_color_manual(values = mal_colors) + scale_x_discrete(labels = mal_labels) +
  theme_bw() + labs(x = NULL, y = "morans_I") + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 18), axis.title = element_text(size = 20))
p16b <- ggplot(morans_compartment, aes(Malignancy, morans_I, color = Malignancy)) +
  geom_boxplot() + stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", label.y = 0.7, step.increase = c(0, 0.05, 0.07), tip.length = 0.01, size = 5) +
  facet_wrap(~cell_type, nrow = 2) + scale_color_manual(values = mal_colors) + scale_x_discrete(labels = mal_labels) +
  theme_bw() + labs(x = NULL, y = "morans_I") + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, size = 15), axis.title = element_text(size = 18), strip.text = element_text(size = 15))

pdf(file.path(plot_dir, "Supplementary_Fig16a.pdf"), width = 4, height = 4); print(p16a); dev.off()
pdf(file.path(plot_dir, "Supplementary_Fig16b.pdf"), width = 7, height = 6); print(p16b); dev.off()

message("Supplementary Figures 14 and 16 complete: ", dirname(out_dir))
