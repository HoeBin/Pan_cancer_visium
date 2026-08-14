# 분석 목적
#   - Visium score 및 LR-Jaccard 중간 결과에서 Figure 4a–d를 재현한다.
#
# 분석 흐름
#   1. compartment annotation과 Visium signature score를 결합한다.
#   2. epithelial·mesenchymal·pEMT score를 sample별로 요약한다.
#   3. Boundary에서 pEMT와 subtype abundance의 Spearman correlation을 계산한다.
#   4. HSPA6+ macrophage–tCAF LR Jaccard 결과를 집계한다.
#
# 주요 출력
#   - Fig4a_a(Epithelial), Fig4a_b(Mesenchymal), Fig4b–d PDF
#   - score, correlation, LR ranking 및 통계 검정 CSV
#
# 출력 위치
#   - Figure/Output/Figure4/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(dplyr); library(tidyr); library(ggplot2); library(ggpubr); library(ggrepel); library(patchwork)})
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
out_dir <- file.path(root_dir, "Output/Figure4/Tables"); plot_dir <- file.path(root_dir, "Output/Figure4/Plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
mal_colors <- c(Mal = "#FF0000", Bdy = "#008B00", Normal = "#4169E1"); mal_labels <- c(Mal = "Malignant", Bdy = "Boundary", Normal = "Normal")
comparisons <- list(c("Normal", "Bdy"), c("Bdy", "Mal"), c("Normal", "Mal"))
celltype_colors <- c(Epithelial = "#FF0000", Tcell = "#0000FF", Bcell = "#00FF00", Myeloid = "#9932CC", Endothelial = "#FF00FF", Fibroblast = "#FFA500", Mural = "#00EEEE")

# "Load score intermediates" --------------------------------------------------
obs <- load_visium_intermediate(include_subtypes = TRUE)
score_file <- file.path(PROCESSED_DATA_DIR, "Visium_Scoring_251214.csv")
scores <- fread(score_file); setnames(scores, names(scores)[1], "cellid")
dat <- merge(obs, scores, by = "cellid", all = FALSE)
score_names <- c("Epithelial", "Mesenchymal", "pEMT_Itaiyanai")
sample_scores <- as.data.frame(dat) %>% group_by(cancer_type, sample_id, Malignancy) %>% summarise(across(all_of(score_names), ~ median(.x, na.rm = TRUE)), .groups = "drop")
write.csv(sample_scores, file.path(out_dir, "Fig4a_4b_sample_signature_scores.csv"), row.names = FALSE)

score_tests <- bind_rows(lapply(score_names, function(score) {
  x <- sample_scores[, c("Malignancy", score)]; names(x)[2] <- "value"
  compare_means(value ~ Malignancy, data = x, method = "wilcox.test", comparisons = comparisons) %>% mutate(signature = score)
}))
write.csv(score_tests, file.path(out_dir, "Fig4a_4b_wilcoxon.csv"), row.names = FALSE)

plot_score <- function(score, title) {
  x <- sample_scores[, c("Malignancy", score)]; names(x)[2] <- "value"
  ggplot(x, aes(Malignancy, value, color = Malignancy)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    scale_color_manual(values = mal_colors) +
    scale_x_discrete(labels = mal_labels) +
    stat_compare_means(method = "wilcox.test", comparisons = comparisons,
                       label = "p.signif", step.increase = c(0, 0.01, 0.03),
                       tip.length = 0.01,
                       label.y = quantile(x$value, 0.75, na.rm = TRUE), size = 6) +
    ylim(c(quantile(x$value, 0.05, na.rm = TRUE),
           quantile(x$value, 0.95, na.rm = TRUE))) +
    theme_bw() +
    labs(title = paste0(title, " (median by sample)"), x = NULL,
         y = "Score (median)") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 15, face = "plain", color = "black"),
          axis.title = element_text(size = 20, color = "black"),
          axis.text = element_text(size = 15, color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1))
}
p_epi <- plot_score("Epithelial", "Epithelial")
p_mes <- plot_score("Mesenchymal", "Mesenchymal")
pdf(file.path(plot_dir, "Fig4a_a.pdf"), width = 3, height = 4)
print(p_epi)
dev.off()
pdf(file.path(plot_dir, "Fig4a_b.pdf"), width = 3, height = 4)
print(p_mes)
dev.off()

sample_pemt <- sample_scores[, c("Malignancy", "pEMT_Itaiyanai")]
names(sample_pemt)[2] <- "value"
p_pemt <- ggplot(sample_pemt, aes(Malignancy, value, color = Malignancy)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  scale_color_manual(values = mal_colors) +
  scale_x_discrete(labels = mal_labels) +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons,
                     label = "p.signif", step.increase = c(0, 0.01, 0.03),
                     tip.length = 0.01,
                     label.y = quantile(sample_pemt$value, 0.75, na.rm = TRUE), size = 6) +
  ylim(c(quantile(sample_pemt$value, 0.05, na.rm = TRUE),
         quantile(sample_pemt$value, 0.95, na.rm = TRUE))) +
  theme_bw() +
  labs(title = "pEMT_Itaiyanai (median by sample)", x = NULL, y = "Score (median)") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15, face = "plain", color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))
pdf(file.path(plot_dir, "Fig4b.pdf"), width = 3, height = 4)
print(p_pemt)
dev.off()

# "Figure 4c: pEMT-abundance correlation" ------------------------------------
metadata <- c("cellid", "Organ", "cancer_type", "sample_id", "Data_GEO", "Sample_GEO", unname(global_celltype_columns()), "Malignancy", "cnv_score", names(scores))
subtype_cols <- setdiff(names(dat), metadata); subtype_cols <- subtype_cols[vapply(dat[, ..subtype_cols], is.numeric, logical(1))]
bdy <- dat[Malignancy == "Bdy"]
subtype_groups <- list(
  Epithelial = paste0("Epi", 0:16),
  Tcell = c("Trm_cytotoxic_T", "Tem/Trm_cytotoxic_T", "CD16-_NK", "exT", "Tem/Effector_helper_T", "regT", "Tem/Temra_cytotoxic_T", "ILC3", "Follicular_helper_T", "Tcm/Naive_helper_T", "FOXP3-_regT", "MAIT", "gamma-delta_T", "Type_17_helper_T", "NK", "Tcm/Naive_cytotoxic_T", "Regulatory_T", "cycling_exT", "CD16+_NK", "ISG15+_T", "cycling_regT", "Type_1_helper_T"),
  Bcell = c("B", "Naive_B", "Memory_B", "Germinal_center_B", "Plasma", "Proliferative_germinal_center_B", "Age-associated_B", "Plasmablasts"),
  Myeloid = c("CXCL3+_Macro", "CD14+_Mono", "CD1C+_cDC", "FOLR2+_Macro", "HSPA6+_Macro", "SPP1+_cycMacro", "CD16+_Mono", "Macro", "LILRA4+_pDC", "CLEC9A+_cDC", "TNFSF10+_Macro", "SPP1+_CXCL3+_Macro", "SPP1+_Macro", "SPP1+_MT1+_Macro", "LAMP3+_cDC", "Mast"),
  Endothelial = c("Venous_EC", "Capillary_EC", "Arterial_EC", "Fibrosis_PGF+_Tip_EC", "Venous_iEC", "PGF+_Tip_EC", "Lymphatic_EC", "HMGB2+_cycEC", "Tip_EC", "SPRY1+_EC", "CCL2+_iEC", "HMOX1+_iEC", "CD14+_circEC", "TMEM100+_EC"),
  Fibroblast = c("CXCL14+_mCAF", "PI16+_iCAF", "iCAF", "Normal_Fibroblast", "apCAF", "vCAF", "IL6+_iCAF", "pnCAF", "mCAF", "ISG15+_mCAF", "tCAF", "HSP+_tCAF", "cyc_mCAF", "myoCAF"),
  Mural = c("SMC", "Pericyte")
)
subtype_to_global <- setNames(rep(names(subtype_groups), lengths(subtype_groups)), unlist(subtype_groups, use.names = FALSE))
celltype_map <- fread(file.path(PROJECT_DIR, "Cell2Location_SubNum4/Data_Info/celltype_mapping.csv"))
rename_key <- setNames(celltype_map$after_celltype, celltype_map$before_celltype)
cor_table <- bind_rows(lapply(subtype_cols, function(column) {
  keep <- complete.cases(bdy[[column]], bdy$pEMT_Itaiyanai)
  test <- suppressWarnings(cor.test(bdy[[column]][keep], bdy$pEMT_Itaiyanai[keep], method = "spearman", exact = FALSE))
  display_name <- if (column %in% names(rename_key)) unname(rename_key[column]) else column
  data.frame(celltype = display_name, global_celltype = unname(subtype_to_global[column]), spearman_rho = unname(test$estimate), p_value = test$p.value, n = sum(keep))
})) %>% filter(!is.na(global_celltype)) %>% arrange(desc(spearman_rho)) %>% mutate(rank = row_number(), global_celltype = factor(global_celltype, levels = names(celltype_colors)))
write.csv(cor_table, file.path(out_dir, "Fig4c_pEMT_abundance_correlations.csv"), row.names = FALSE)
labelled_subtypes <- c(
  "cyc_mCAF", "ISG15+_mCAF", "mCAF", "tCAF", "SPP1+_cycMacro",
  "SPP1+_MT1+_Macro", "SPP1+_CXCL3+_Macro", "Epi_Interferon",
  "SPP1+_Macro", "HSPA6+_Macro"
)
top_cor <- cor_table %>%
  filter(celltype %in% labelled_subtypes) %>%
  mutate(celltype = factor(celltype, levels = labelled_subtypes)) %>%
  arrange(celltype) %>%
  mutate(plot_label = as.character(celltype))
p4c <- ggplot(cor_table, aes(rank, spearman_rho)) + geom_point(aes(color = global_celltype), size = 2) + geom_point(data = top_cor, size = 3.5, color = "grey20") + geom_point(data = top_cor, aes(color = global_celltype), size = 2) + geom_text_repel(data = top_cor, aes(label = plot_label), size = 6, box.padding = 0.8, point.padding = 2, max.overlaps = Inf, force = 5, force_pull = 0.1, min.segment.length = 0, segment.curvature = 0.1, segment.alpha = 0.3) + scale_color_manual(values = celltype_colors, guide = guide_legend(override.aes = list(size = 6)), labels = c(Epithelial = "Epithelial", Tcell = "T cell", Bcell = "B cell", Myeloid = "Myeloid", Endothelial = "Endothelial", Fibroblast = "Fibroblast", Mural = "Mural")) + theme_bw() + labs(title = "pEMT_Itaiyanai", x = NULL, y = "Correlation", color = "celltype") + theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size = 20, face = "plain"), axis.title = element_text(size = 15), axis.text = element_text(color = "black", size = 15), axis.text.x = element_blank(), axis.ticks.x = element_blank())
pdf(file.path(plot_dir, "Fig4c.pdf"), width = 9, height = 5); print(p4c); dev.off()

# "Figure 4d: HSPA6+ macrophage-tCAF LR pairs" -------------------------------
lr_root <- file.path(LR_DIR, "LR_Sub_All_CT_Jaccard/Output/All_Subtypes_bdy")
lr_base <- "GroupPair_LR_Jaccard_each_sample_Myeloid@Fibroblast"
lr_data <- bind_rows(lapply(c("Mal", "Bdy", "Normal"), function(region) {
  path <- file.path(lr_root, paste0(lr_base, "_", region, ".csv"))
  cmd <- sprintf("LC_ALL=C grep -E %s %s", shQuote("^(HSPA6\\+_Macro@tCAF|tCAF@HSPA6\\+_Macro),"), shQuote(path))
  x <- fread(cmd = cmd, header = FALSE, showProgress = FALSE)
  setnames(x, c("group_pair", "lr_pair", "cancer_type", "sample_id", "LocationNew", "jaccard")); x$Region <- region; x
}))
lr_data <- lr_data[(grepl("HSPA6[+]_Macro", group_pair) & grepl("tCAF", group_pair)) & sample_id %in% unique(obs$sample_id)]
pair_per_sample <- as.data.frame(lr_data) %>% count(sample_id, lr_pair, name = "n")
common_pairs <- pair_per_sample %>% group_by(lr_pair) %>% summarise(sample_count = n_distinct(sample_id), .groups = "drop") %>% filter(sample_count == n_distinct(lr_data$sample_id)) %>% pull(lr_pair)
lr_common <- as.data.frame(lr_data) %>% filter(lr_pair %in% common_pairs)
lr_summary <- lr_common %>% group_by(lr_pair, Region) %>% summarise(median_jaccard = median(jaccard, na.rm = TRUE), mean_jaccard = mean(jaccard, na.rm = TRUE), n = n(), .groups = "drop") %>% mutate(plot_value = mean_jaccard, Region = factor(Region, levels = c("Mal", "Bdy", "Normal")))
lr_rank <- lr_summary %>% select(lr_pair, Region, plot_value) %>% pivot_wider(names_from = Region, values_from = plot_value) %>% mutate(mean_MN = (Mal + Normal) / 2, max_MN = pmax(Mal, Normal), boundary_gain = Bdy - max_MN) %>% arrange(desc(boundary_gain)) %>% mutate(rank = row_number())
write.csv(lr_data, file.path(out_dir, "Fig4d_LR_plot_data.csv"), row.names = FALSE); write.csv(lr_summary, file.path(out_dir, "Fig4d_LR_summary_by_region.csv"), row.names = FALSE); write.csv(lr_rank, file.path(out_dir, "Fig4d_boundary_LR_rank.csv"), row.names = FALSE)
paper_bottom_pairs <- c("ALOX5AP-ALOX5", "S100A9-ITGB2")
lr_bdy_rank <- lr_summary %>%
  filter(Region == "Bdy", median_jaccard > 0) %>%
  arrange(desc(median_jaccard)) %>%
  mutate(rank = row_number())
label_pairs <- unique(c(slice_head(lr_bdy_rank, n = 8)$lr_pair, paper_bottom_pairs))
lr_bdy_labels <- filter(lr_bdy_rank, lr_pair %in% label_pairs)
write.csv(lr_bdy_rank, file.path(out_dir, "Fig4d_paper_boundary_median_rank.csv"), row.names = FALSE)
p4d <- ggplot(lr_bdy_rank, aes(rank, median_jaccard)) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(color = "black", size = 1.7) +
  geom_point(data = lr_bdy_labels, color = "red", size = 3) +
  geom_text_repel(data = lr_bdy_labels, aes(label = lr_pair),
                  nudge_x = 90, direction = "y", hjust = 0,
                  size = 3.5, box.padding = 0.5, point.padding = 0.25,
                  min.segment.length = 0, segment.color = "grey30",
                  force = 5, force_pull = 0, max.overlaps = Inf,
                  xlim = c(70, 175)) +
  scale_x_continuous(limits = c(0, 220), expand = expansion(mult = c(0.02, 0.02))) +
  theme_bw() +
  labs(title = "Boundary co-localization between\nHSPA6+ Macro–tCAF and LR pairs",
       x = "Rank", y = "Median Co-localization score") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "plain"),
        axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(5.5, 30, 5.5, 5.5))
pdf(file.path(plot_dir, "Fig4d.pdf"), width = 6, height = 3.6)
print(p4d)
dev.off()
message("Figure 4 complete: ", dirname(out_dir))
