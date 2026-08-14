# 분석 목적
#   - correlation, Jaccard 및 network centrality 중간 결과에서 Figure 3 패널을 재현한다.
#
# 분석 흐름
#   1. Global correlation과 암종별 co-localization summary를 시각화한다.
#   2. 3-region subtype Jaccard 결과로 compartment 비교를 수행한다.
#   3. Boundary co-enrichment network의 edge/node table을 정리한다.
#   4. Cytoscape betweenness 결과로 global median 및 Boundary rank를 계산한다.
#   5. 주요 B-cell/macrophage subtype pair를 선택해 boxplot을 생성한다.
#
# 주요 출력
#   - Fig3a,b,c,e,f,g,h PDF
#   - Fig3d network edge/node table (Cytoscape 원본 다이어그램이라 PDF 재현 대상 아님)
#   - 패널별 plot data, summary, Wilcoxon CSV
#
# 출력 위치
#   - Figure/Output/Figure3/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(dplyr); library(tidyr); library(ggplot2); library(ggpubr); library(patchwork); library(ComplexHeatmap); library(circlize)})
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
out_dir <- file.path(root_dir, "Output/Figure3/Tables"); plot_dir <- file.path(root_dir, "Output/Figure3/Plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
project <- PROJECT_DIR
colors <- c(Epithelial = "red", `T cell` = "blue", `B cell` = "green", Myeloid = "purple", Endothelial = "magenta", Fibroblast = "orange", Mural = "cyan")
mal_colors <- c(Mal = "#FF0000", Bdy = "#008B00", Normal = "#4169E1"); mal_labels <- c(Mal = "Malignant", Bdy = "Boundary", Normal = "Normal")
comparisons <- list(c("Normal", "Bdy"), c("Bdy", "Mal"), c("Normal", "Mal"))

# "Figure 3a: global correlation heatmap" ------------------------------------
obs_global <- load_visium_intermediate(include_subtypes = FALSE)
global_cols <- unname(global_celltype_columns())
mat <- cor(as.matrix(obs_global[, ..global_cols]), use = "pairwise.complete.obs")
display_names <- names(global_celltype_columns())
rownames(mat) <- display_names
colnames(mat) <- display_names
cor_summary <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
names(cor_summary) <- c("celltype_1", "celltype_2", "correlation")
write.csv(cor_summary, file.path(out_dir, "Fig3a_global_correlation.csv"), row.names = FALSE)
ht3a <- Heatmap(
  mat,
  column_title = "All Cancer Types",
  show_row_dend = TRUE,
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  row_names_gp = gpar(col = "black", fontsize = 15),
  column_names_gp = gpar(col = "black", fontsize = 15),
  column_title_gp = gpar(fontsize = 20, fontface = "bold")
)
pdf(file.path(plot_dir, "Fig3a.pdf"), width = 4, height = 4)
draw(ht3a, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
dev.off()

# "Figure 3b: cancer-wise global pair co-localization" ------------------------
j3b <- fread(file.path(project, "Cell2Location_SubNum4/Supple_Table/Fig3b_jaccard.csv")); j3b <- j3b[cancer_type != "ESCA"]
write.csv(j3b, file.path(out_dir, "Fig3b_key_results_jaccard_by_cancer.csv"), row.names = FALSE)
pair_order <- j3b %>% group_by(Feat_pair) %>% summarise(x = median(median, na.rm = TRUE), .groups = "drop") %>% arrange(desc(x)) %>% pull(Feat_pair)
p3b <- ggplot(j3b, aes(factor(Feat_pair, levels = pair_order), median)) + geom_boxplot() + geom_point(aes(color = cancer_type, shape = cancer_type), position = position_jitter(width = 0.12), size = 2) + theme_bw() + labs(x = NULL, y = "Median Co-localization Score") + theme(axis.text.x = element_blank(), legend.position = "right")
pdf(file.path(plot_dir, "Fig3b.pdf"), width = 10, height = 5); print(p3b); dev.off()

# "Figure 3c and subtype panels" ---------------------------------------------
sub_raw <- fread(file.path(STOPOVER_DIR, "subtype/subtype_3region.csv"))
non_esca_samples <- unique(obs_global$sample_id)
sub_raw[, Sample_clean := sub("^sp_concat_DEG_250618__", "", Sample)]
sub_raw <- sub_raw[Sample_clean %in% non_esca_samples]
sub_raw[, Sample_clean := NULL]
sub_raw$Region <- factor(sub_raw$Region, levels = c("Mal", "Bdy", "Normal"))
sub <- copy(sub_raw)
celltype_map <- fread(file.path(project, "Cell2Location_SubNum4/Data_Info/celltype_mapping.csv"))
rename_key <- setNames(celltype_map$after_celltype, celltype_map$before_celltype)
sub[, Feat_1 := fifelse(Feat_1 %chin% names(rename_key), rename_key[Feat_1], Feat_1)]
sub[, Feat_2 := fifelse(Feat_2 %chin% names(rename_key), rename_key[Feat_2], Feat_2)]
sub[, Feat_pair := paste0(Feat_1, "@", Feat_2)]
sub_summary <- sub[, .(
  median = median(J_comp, na.rm = TRUE),
  mean = mean(J_comp, na.rm = TRUE),
  sd = sd(J_comp, na.rm = TRUE)
), by = .(Region, Feat_pair, Feat1 = Feat_1, Feat2 = Feat_2)]
sub_summary$Region <- factor(sub_summary$Region, levels = c("Mal", "Bdy", "Normal"))
write.csv(sub_summary, file.path(out_dir, "Fig3c_subtype_jaccard_plot_data.csv"), row.names = FALSE)
tests3c <- compare_means(median ~ Region, data = sub_summary, method = "wilcox.test", comparisons = comparisons)
write.csv(tests3c, file.path(out_dir, "Fig3c_wilcoxon.csv"), row.names = FALSE)
p3c <- ggplot(sub_summary, aes(Region, median)) +
  geom_jitter(size = 0.5, fill = "black", alpha = 0.1) +
  geom_boxplot(aes(color = Region), alpha = 0.8, outlier.shape = NA) +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons,
                     label = "p.signif", step.increase = c(0.05, 0.05, 0.08),
                     label.y = 0.75, size = 6) +
  scale_x_discrete(labels = mal_labels) +
  scale_color_manual(values = mal_colors, labels = mal_labels, name = "Malignancy") +
  theme_bw() +
  labs(x = NULL, y = "Jaccard Index") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold", color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 20, face = "plain"),
        legend.text = element_text(size = 15, color = "black"))
pdf(file.path(plot_dir, "Fig3c.pdf"), width = 6, height = 4); print(p3c); dev.off()

plot_target <- function(target, panel) {
  dat <- sub[Feat_1 == target | Feat_2 == target]
  dat$partner <- ifelse(dat$Feat_1 == target, dat$Feat_2, dat$Feat_1)
  top <- as.data.frame(dat) %>% filter(Region == "Bdy") %>% group_by(partner) %>% summarise(x = median(J_comp, na.rm = TRUE), .groups = "drop") %>% arrange(desc(x)) %>% slice_head(n = 10) %>% pull(partner)
  dat <- dat[partner %in% top]; dat$partner <- factor(dat$partner, levels = top)
  write.csv(dat, file.path(out_dir, paste0(panel, "_plot_data.csv")), row.names = FALSE)
  test <- compare_means(J_comp ~ Region, data = dat, group.by = "partner", method = "wilcox.test", comparisons = comparisons)
  write.csv(test, file.path(out_dir, paste0(panel, "_wilcoxon.csv")), row.names = FALSE)
  p <- ggplot(dat, aes(Region, J_comp, fill = Region)) +
    geom_boxplot(width = 1) +
    facet_wrap(~partner, nrow = 1, strip.position = "bottom") +
    stat_compare_means(method = "wilcox.test", comparisons = list(c("Bdy", "Mal"), c("Normal", "Bdy"), c("Normal", "Mal")), label = "p.signif", step.increase = c(0.06, 0.06, 0.1), label.y = 1, size = 6) +
    ylim(0, 1.3) + theme_bw() +
    scale_x_discrete(labels = mal_labels) + scale_fill_manual(values = mal_colors, breaks = c("Mal", "Bdy", "Normal"), labels = mal_labels, name = "Malignancy") +
    labs(x = NULL, y = "Median Jaccard Index") +
    theme(strip.background = element_blank(), strip.text.x = element_text(size = 20, hjust = 1, vjust = 1.05, angle = 45, color = "black", margin = margin(t = 10)), strip.placement = "outside", panel.spacing = unit(4, "mm"), strip.clip = "off", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 15, color = "black"), axis.title = element_text(size = 20, color = "black"), legend.position = "bottom", legend.direction = "horizontal", plot.margin = margin(20, 20, 20, 20, unit = "mm"))
  pdf(file.path(plot_dir, paste0(panel, ".pdf")), width = 11.5, height = 9); print(p); dev.off()
}
plot_target("Cycling_Bgc", "Fig3g")
plot_target("HSPA6+_Macro", "Fig3h")

# "Figure 3d: boundary co-enrichment network tables" --------------------------
# Cytoscape에서 수기로 배치한 다이어그램이라 plot 자체는 재현 대상이 아니며,
# 입력으로 쓰인 edge(Jaccard co-localization score > 0.5 threshold)/node table만 정리한다.
node_dir <- file.path(project, "Cell2Location_SubNum4/Network/Stopover2")
edge3d <- fread(file.path(node_dir, "Cor_cluster_Bdy_cut0.5.csv"))
edge3d[, Var1 := fifelse(Var1 %chin% names(rename_key), rename_key[Var1], Var1)]
edge3d[, Var2 := fifelse(Var2 %chin% names(rename_key), rename_key[Var2], Var2)]
edge3d[, Var3 := paste(Var1, "(interacts with)", Var2)]
write.csv(edge3d, file.path(out_dir, "Fig3d_boundary_network_edges.csv"), row.names = FALSE)
node3d <- fread(file.path(node_dir, "Cor_cluster_Bdy_cut0.5.csv_1 default node.csv"))
node3d[, name := fifelse(name %chin% names(rename_key), rename_key[name], name)]
node3d[, `shared name` := name]
write.csv(node3d, file.path(out_dir, "Fig3d_boundary_network_nodes.csv"), row.names = FALSE)

# "Figure 3e-f: betweenness centrality" --------------------------------------
bc <- fread(file.path(project, "Cell2Location_SubNum4/Supple_Table/Fig3e_BC.csv"))
present_names_raw <- unique(unlist(lapply(c("Mal", "Bdy", "Normal"), function(region) {
  fread(file.path(node_dir, paste0("Cor_cluster_", region, "_cut0.5.csv_1 default node.csv")), select = "name")$name
})))
present_names <- ifelse(present_names_raw %in% names(rename_key), rename_key[present_names_raw], present_names_raw)
bc_for_summary <- bc %>% filter(name %in% present_names)
bc_long <- bc_for_summary %>% pivot_longer(c(Mal, Bdy, Normal), names_to = "Region", values_to = "BetweennessCentrality")
bc_summary <- bc_long %>% group_by(global_cell_type, Region) %>% summarise(median_betweenness = median(BetweennessCentrality, na.rm = TRUE), mean_betweenness = mean(BetweennessCentrality, na.rm = TRUE), n_subtypes = n(), .groups = "drop")
write.csv(bc, file.path(out_dir, "Fig3e_3f_betweenness_by_subtype.csv"), row.names = FALSE); write.csv(bc_summary, file.path(out_dir, "Fig3e_key_results_median_betweenness.csv"), row.names = FALSE)
global_order <- c("Epithelial", "T cell", "B cell", "Myeloid", "Endothelial", "Fibroblast", "Mural")
bc_summary <- bc_summary %>%
  mutate(global_cell_type = factor(global_cell_type, levels = global_order),
         Region = factor(Region, levels = c("Mal", "Bdy", "Normal")))
showtext::showtext_auto()
p3e_all <- ggplot(bc_summary, aes(Region, median_betweenness, fill = Region)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~global_cell_type, scales = "free_y") + theme_bw() +
  scale_fill_manual(values = mal_colors) + scale_x_discrete(labels = mal_labels) +
  labs(x = NULL, y = "Median Betweenness Centrality") +
  theme(legend.position = "none", legend.title = element_text(size = 30),
        legend.text = element_text(size = 20), legend.key.size = unit(1.5, "cm"),
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.title = element_text(hjust = 0.5, size = 20),
        axis.text = element_text(color = "black", size = 15),
        axis.text.x = element_text(size = 20, angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        strip.text = element_text(color = "black", size = 20)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.0001))
p3e_bdy <- filter(bc_summary, Region == "Bdy") %>%
  ggplot(aes(global_cell_type, median_betweenness, fill = global_cell_type)) +
  geom_col(position = "dodge", color = "black") + theme_bw() +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = NULL, title = "Median Betweenness Centrality") +
  theme(legend.position = "none", legend.title = element_text(size = 30),
        legend.text = element_text(size = 20), legend.key.size = unit(1.5, "cm"),
        plot.title = element_text(hjust = 0.5, size = 20, face = "plain"),
        axis.title = element_text(hjust = 0.5, size = 20),
        axis.text = element_text(color = "black", size = 15),
        axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(color = "black", size = 20)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001))
pdf(file.path(plot_dir, "Fig3e.pdf"), width = 7, height = 7)
print(p3e_all)
print(p3e_bdy)
dev.off()
bc_rank <- bc %>% arrange(desc(Bdy)) %>% mutate(rank = row_number())
write.csv(bc_rank, file.path(out_dir, "Fig3f_boundary_betweenness_rank.csv"), row.names = FALSE)
p3f <- ggplot(bc_rank, aes(rank, Bdy, color = global_cell_type)) + geom_point(size = 2) + geom_text(data = slice_head(bc_rank, n = 12), aes(label = name), hjust = 0, nudge_x = 2, size = 3) + scale_color_manual(values = colors) + theme_bw() + labs(title = "Boundary", x = "Rank", y = "Betweenness Centrality")
pdf(file.path(plot_dir, "Fig3f.pdf"), width = 6, height = 4); print(p3f); dev.off()
message("Figure 3 complete: ", dirname(out_dir))
