# 분석 목적
#   - Visium 중간 결과에서 Supplementary Figure 8의 CNV 관련 패널을 재현한다.
#
# 분석 흐름
#   1. 15개 Obs.csv와 malignancy/CNV 결과를 결합한다.
#   2. 영역별 CNV 및 spot 수를 요약하고 Wilcoxon 검정을 수행한다.
#   3. Global/subtype proportion과 CNV score의 상관계수를 계산한다.
#   4. 패널 PDF와 Google Sheets용 CSV를 저장한다.
#
# 주요 출력
#   - Supplementary_Fig8a–d PDF
#   - plot data, correlation, Wilcoxon 결과 CSV
#
# 출력 위치
#   - Figure/Output/Supplementary_Figure8/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  library(scales)
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

out_dir  <- file.path(root_dir, "Output/Supplementary_Figure8/Tables")
plot_dir <- file.path(root_dir, "Output/Supplementary_Figure8/Plots")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# "Load plot-ready intermediate data" ----------------------------------------
obs <- load_visium_intermediate(include_subtypes = TRUE)
global_map <- global_celltype_columns()
mal_colors <- c(Mal = "#FF0000", Bdy = "#008B00", Normal = "#4169E1")
mal_labels <- c(Mal = "Malignant", Bdy = "Boundary", Normal = "Normal")
comparisons <- list(c("Mal", "Bdy"), c("Bdy", "Normal"), c("Mal", "Normal"))

# "Panel 8a: CNV score" -------------------------------------------------------
cnv_slide <- obs[, .(
  mean_cnv = mean(cnv_score, na.rm = TRUE),
  median_cnv = median(cnv_score, na.rm = TRUE),
  n_spots = .N
), by = .(cancer_type, sample_id, Malignancy)]

write.csv(obs[, .(cellid, cancer_type, sample_id, Malignancy, cnv_score)],
          file.path(out_dir, "Supplementary_Fig8a_spot_level_cnv.csv"), row.names = FALSE)
write.csv(cnv_slide, file.path(out_dir, "Supplementary_Fig8a_slide_level_cnv.csv"), row.names = FALSE)

cnv_tests <- compare_means(cnv_score ~ Malignancy, data = as.data.frame(obs),
                           method = "wilcox.test", comparisons = comparisons)
cnv_median_tests <- compare_means(median_cnv ~ Malignancy, data = as.data.frame(cnv_slide),
                                  method = "wilcox.test", comparisons = comparisons)
write.csv(cnv_tests, file.path(out_dir, "Supplementary_Fig8a_cnv_wilcoxon.csv"), row.names = FALSE)
write.csv(cnv_median_tests, file.path(out_dir, "Supplementary_Fig8a_median_cnv_wilcoxon.csv"), row.names = FALSE)

p8a_spot <- ggplot(obs, aes(Malignancy, cnv_score, color = Malignancy)) +
  geom_boxplot() +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                     label.y = c(6500, 6600, 7600), label = "p.signif") +
  scale_color_manual(values = mal_colors) +
  scale_x_discrete(labels = mal_labels) +
  scale_y_continuous(labels = comma_format(accuracy = 1)) +
  theme_cowplot() + labs(title = "CNV score", x = NULL, y = "CNV") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 13))

p8a_slide <- ggplot(cnv_slide, aes(Malignancy, median_cnv, color = Malignancy)) +
  geom_jitter(color = "black", alpha = 0.3, size = 0.5) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                     label.y = c(6500, 6600, 7600), label = "p.signif") +
  scale_color_manual(values = mal_colors) +
  scale_x_discrete(labels = mal_labels) +
  scale_y_continuous(labels = comma_format(accuracy = 1)) +
  theme_cowplot() + labs(title = "CNV score", x = NULL, y = "CNV") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 13))

pdf(file.path(plot_dir, "Supplementary_Fig8a.pdf"), width = 4, height = 3)
print(p8a_spot)
print(p8a_slide)
dev.off()

# "Panel 8b: number of spots" -------------------------------------------------
spot_counts <- obs[, .(n = .N), by = .(cancer_type, sample_id, Malignancy)]
spot_tests <- compare_means(n ~ Malignancy, data = as.data.frame(spot_counts),
                            method = "wilcox.test", comparisons = comparisons)
write.csv(spot_counts, file.path(out_dir, "Supplementary_Fig8b_spot_counts.csv"), row.names = FALSE)
write.csv(spot_tests, file.path(out_dir, "Supplementary_Fig8b_spot_count_wilcoxon.csv"), row.names = FALSE)

p8b <- ggplot(spot_counts, aes(Malignancy, n, color = Malignancy)) +
  geom_jitter(color = "black", alpha = 0.3, size = 0.5) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                     label.y = c(3000, 3300, 3900), label = "p.signif") +
  scale_color_manual(values = mal_colors) +
  scale_x_discrete(labels = mal_labels) +
  scale_y_continuous(labels = comma_format(accuracy = 1)) +
  theme_cowplot() + labs(title = "Num of Spot", x = NULL, y = "Count") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

pdf(file.path(plot_dir, "Supplementary_Fig8b.pdf"), width = 4, height = 3)
print(p8b)
dev.off()

# "Panels 8c-d: CNV correlations" --------------------------------------------
cor_one <- function(column, label) {
  test <- suppressWarnings(cor.test(obs$cnv_score, obs[[column]], method = "pearson"))
  data.frame(
    celltype = label,
    source_column = column,
    correlation = unname(test$estimate),
    p_value = test$p.value,
    n = sum(stats::complete.cases(obs$cnv_score, obs[[column]]))
  )
}

global_cor <- bind_rows(Map(cor_one, unname(global_map), names(global_map)))
metadata_cols <- c(
  "cellid", "Organ", "cancer_type", "sample_id", "Data_GEO", "Sample_GEO",
  unname(global_map), "Malignancy", "cnv_score"
)
subtype_cols <- setdiff(names(obs), metadata_cols)
subtype_cols <- subtype_cols[vapply(obs[, ..subtype_cols], is.numeric, logical(1))]

celltype_map <- fread(file.path(PROJECT_DIR, "Cell2Location_SubNum4/Data_Info/celltype_mapping.csv"))
rename_key <- setNames(celltype_map$after_celltype, celltype_map$before_celltype)
subtype_cor <- bind_rows(lapply(subtype_cols, function(x) {
  display_name <- if (x %in% names(rename_key)) unname(rename_key[x]) else x
  cor_one(x, display_name)
}))

write.csv(global_cor, file.path(out_dir, "Supplementary_Fig8c_global_cnv_correlations.csv"), row.names = FALSE)
write.csv(subtype_cor, file.path(out_dir, "Supplementary_Fig8d_subtype_cnv_correlations.csv"), row.names = FALSE)

global_colors <- c(
  Epithelial = "red", Tcell = "blue", Bcell = "green",
  Myeloid = "purple", Endothelial = "magenta",
  Fibroblast = "orange", Mural = "cyan"
)
p8c <- ggplot(global_cor, aes(reorder(celltype, -correlation), correlation, fill = celltype)) +
  geom_col(color = "black") + geom_hline(yintercept = 0) +
  scale_fill_manual(values = global_colors) + theme_cowplot() +
  labs(x = NULL, y = "Correlation", title = "Correlation with global proportion and CNV score") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        plot.title = element_text(hjust = 0.5, size = 13))

major_category_map <- setNames(celltype_map$major_category, celltype_map$before_celltype)
major_to_parent <- c(Epithelial = "Epithelial", T_NK = "Tcell", B_cell = "Bcell", Myeloid = "Myeloid", Endothelial = "Endothelial", Fibroblast = "Fibroblast", Mural = "Mural", Unknown = "Fibroblast")
subtype_cor$parent <- unname(major_to_parent[major_category_map[subtype_cor$source_column]])
p8d <- ggplot(subtype_cor, aes(reorder(celltype, -correlation), correlation, fill = parent)) +
  geom_col(color = "black") + geom_hline(yintercept = 0) +
  scale_fill_manual(values = global_colors) + theme_cowplot() +
  labs(x = NULL, y = "Correlation", title = "Correlation with sub proportion and CNV score") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
        plot.title = element_text(hjust = 0.5, size = 13))

pdf(file.path(plot_dir, "Supplementary_Fig8c.pdf"), width = 4, height = 3)
print(p8c)
dev.off()
pdf(file.path(plot_dir, "Supplementary_Fig8d.pdf"), width = 10, height = 4)
print(p8d)
dev.off()

message("Supplementary Figure 8 complete: ", dirname(out_dir))
