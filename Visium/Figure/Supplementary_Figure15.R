# 분석 목적
#   - Stopover Jaccard 중간 결과에서 Supplementary Fig.15를 재현한다.
#
# 분석 흐름
#   1. Global/subtype 3-region Jaccard 결과를 불러온다.
#   2. Global pair는 암종별 median Jaccard를 계산한다.
#   3. Subtype pair는 Epi–Epi·Epi–TME·TME–TME로 분류한다.
#   4. plot data와 Wilcoxon 검정 결과를 CSV로 저장한다.
#
# 주요 출력
#   - Supplementary_Fig15a,b PDF
#   - Global/subtype Jaccard summary 및 Wilcoxon CSV
#
# 출력 위치
#   - Figure/Output/Supplementary_Figure15/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
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
out_dir <- file.path(root_dir, "Output/Supplementary_Figure15/Tables")
plot_dir <- file.path(root_dir, "Output/Supplementary_Figure15/Plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

global_file <- file.path(STOPOVER_DIR, "global/3region_reslt.csv")
subtype_file <- file.path(STOPOVER_DIR, "subtype/subtype_3region.csv")
mal_colors <- c(Mal = "#FF0000", Bdy = "#008B00", Normal = "#4169E1")
mal_labels <- c(Mal = "Malignant", Bdy = "Boundary", Normal = "Normal")
comparisons <- list(c("Normal", "Bdy"), c("Bdy", "Mal"), c("Normal", "Mal"))

# "Sample-to-cancer mapping" --------------------------------------------------
obs <- load_visium_intermediate(include_subtypes = FALSE)
sample_info <- unique(as.data.frame(obs)[, c("sample_id", "cancer_type")])

# "Panel 15a: global pairs" ---------------------------------------------------
global <- fread(global_file)
global$Sample <- sub("^sp_concat_DEG_250618__", "", global$Sample)
global <- merge(global, sample_info, by.x = "Sample", by.y = "sample_id", all = FALSE)
global$Feat_pair <- paste(global$Feat_1, global$Feat_2, sep = "-")
global_summary <- global %>% group_by(Feat_pair, cancer_type, Region) %>%
  summarise(median_J_comp = median(J_comp, na.rm = TRUE), mean_J_comp = mean(J_comp, na.rm = TRUE), n = n(), .groups = "drop")
global_summary$Region <- factor(global_summary$Region, levels = c("Mal", "Bdy", "Normal"))
write.csv(global_summary, file.path(out_dir, "Supplementary_Fig15a_global_pair_jaccard_by_cancer.csv"), row.names = FALSE)

global_tests <- compare_means(median_J_comp ~ Region, data = global_summary, group.by = "Feat_pair", method = "wilcox.test", comparisons = comparisons)
write.csv(global_tests, file.path(out_dir, "Supplementary_Fig15a_wilcoxon_by_pair.csv"), row.names = FALSE)

pair_order <- global_summary %>% filter(Region == "Bdy") %>% group_by(Feat_pair) %>% summarise(x = median(median_J_comp), .groups = "drop") %>% arrange(desc(x)) %>% pull(Feat_pair)
p15a_top <- ggplot(global_summary, aes(factor(Feat_pair, levels = pair_order), median_J_comp, fill = Region)) +
  geom_boxplot(outlier.alpha = 0) + scale_fill_manual(values = mal_colors, labels = mal_labels) +
  theme_bw() + coord_cartesian(ylim = c(0, 1)) + labs(x = NULL, y = "J_comp", title = "J_comp (Global Annotation)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), axis.title = element_text(size = 18), legend.title = element_text(size = 14))

pair_parts <- do.call(rbind, strsplit(pair_order, "-", fixed = TRUE))
membership <- data.frame(Pair = rep(pair_order, each = 2), CellType = as.vector(t(pair_parts)), stringsAsFactors = FALSE)
membership$Pair <- factor(membership$Pair, levels = pair_order)
cell_colors <- c(Epithelial = "red", `T cell` = "blue", `B cell` = "green", Myeloid = "purple", Endothelial = "magenta", Fibroblast = "orange", Mural = "cyan")
p15a_matrix <- ggplot(membership, aes(Pair, CellType, group = Pair)) + geom_line(color = "grey40", linewidth = 2) + geom_point(aes(color = CellType), size = 4) + scale_color_manual(values = cell_colors) + theme_bw() + labs(x = NULL, y = "Cell Type Pair") + theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")
pdf(file.path(plot_dir, "Supplementary_Fig15a.pdf"), width = 13, height = 6)
print(p15a_top / p15a_matrix + plot_layout(heights = c(3, 1)))
dev.off()

# "Panel 15b: subtype pair classes" ------------------------------------------
subtype <- fread(subtype_file)
celltype_map <- fread(file.path(PROJECT_DIR, "Cell2Location_SubNum4/Data_Info/celltype_mapping.csv"))
rename_map <- setNames(celltype_map$after_celltype, celltype_map$before_celltype)
subtype$Feat_1 <- ifelse(subtype$Feat_1 %in% names(rename_map), unname(rename_map[subtype$Feat_1]), subtype$Feat_1)
subtype$Feat_2 <- ifelse(subtype$Feat_2 %in% names(rename_map), unname(rename_map[subtype$Feat_2]), subtype$Feat_2)
subtype$Feat_pair <- paste(subtype$Feat_1, subtype$Feat_2, sep = "@")
subtype_summary <- subtype %>% group_by(Region, Feat_pair) %>% summarise(median_J_comp = median(J_comp, na.rm = TRUE), mean_J_comp = mean(J_comp, na.rm = TRUE), n = n(), .groups = "drop") %>%
  separate(Feat_pair, into = c("Feat1", "Feat2"), sep = "@", remove = FALSE)
subtype_summary$Feat1_type <- ifelse(grepl("^Epi_", subtype_summary$Feat1), "Epi", "TME")
subtype_summary$Feat2_type <- ifelse(grepl("^Epi_", subtype_summary$Feat2), "Epi", "TME")
subtype_summary$Feat_type_pair <- paste(subtype_summary$Feat1_type, subtype_summary$Feat2_type, sep = "-")
subtype_summary$Feat_type_pair[subtype_summary$Feat_type_pair == "TME-Epi"] <- "Epi-TME"
subtype_summary$Region <- factor(subtype_summary$Region, levels = c("Mal", "Bdy", "Normal"))
write.csv(subtype_summary, file.path(out_dir, "Supplementary_Fig15b_subtype_pair_jaccard.csv"), row.names = FALSE)
subtype_tests <- compare_means(median_J_comp ~ Region, data = subtype_summary, group.by = "Feat_type_pair", method = "wilcox.test", comparisons = comparisons)
write.csv(subtype_tests, file.path(out_dir, "Supplementary_Fig15b_wilcoxon_by_pair_class.csv"), row.names = FALSE)

p15b <- ggplot(subtype_summary, aes(Region, median_J_comp, color = Region)) +
  geom_jitter(size = 0.5, alpha = 0.1, position = position_jitter(width = 0.15)) + geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~Feat_type_pair, scales = "free_y") + scale_color_manual(values = mal_colors, labels = mal_labels) +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif", step.increase = c(0.05, 0.05, 0.08), tip.length = 0, size = 5) +
  scale_x_discrete(labels = mal_labels) + theme_bw() + labs(x = NULL, y = "Jaccard Index") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13), axis.title = element_text(size = 18), strip.text = element_text(size = 15), legend.position = "none")
pdf(file.path(plot_dir, "Supplementary_Fig15b.pdf"), width = 8, height = 4)
print(p15b)
dev.off()

message("Supplementary Figure 15 complete: ", dirname(out_dir))
