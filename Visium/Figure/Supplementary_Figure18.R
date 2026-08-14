# 분석 목적
#   - Cytoscape node-table 중간 결과에서 Supplementary Fig.18을 재현한다.
#
# 분석 흐름
#   1. Malignant·Boundary·Normal network node table을 결합한다.
#   2. subtype 명칭과 global cell-type 그룹을 정리한다.
#   3. subtype별 betweenness centrality bar plot을 생성한다.
#   4. 전체 값과 Boundary 증가 순위를 CSV로 저장한다.
#
# 주요 출력
#   - Supplementary_Fig18.pdf
#   - subtype별 centrality 및 Boundary 증가 순위 CSV
#
# 출력 위치
#   - Figure/Output/Supplementary_Figure18/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
root_dir <- dirname(script_path)
source(file.path(root_dir, "common/config.R"))
out_dir <- file.path(root_dir, "Output/Supplementary_Figure18/Tables")
plot_dir <- file.path(root_dir, "Output/Supplementary_Figure18/Plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

network_dir <- file.path(PROJECT_DIR, "Cell2Location_SubNum4/Network/Stopover2")
files <- c(
  Mal = "Cor_cluster_Mal_cut0.5.csv_1 default node.csv",
  Bdy = "Cor_cluster_Bdy_cut0.5.csv_1 default node.csv",
  Normal = "Cor_cluster_Normal_cut0.5.csv_1 default node.csv"
)
celltype_map <- fread(file.path(PROJECT_DIR, "Cell2Location_SubNum4/Data_Info/celltype_mapping.csv"))
rename_map <- setNames(celltype_map$after_celltype, celltype_map$before_celltype)

# "Load network intermediates" ------------------------------------------------
net <- bind_rows(lapply(names(files), function(category) {
  x <- fread(file.path(network_dir, files[[category]]))
  x$Category <- category
  x
}))
net$name <- ifelse(net$name %in% names(rename_map), unname(rename_map[net$name]), net$name)
score <- net %>% select(name, Category, Score = BetweennessCentrality) %>%
  pivot_wider(names_from = Category, values_from = Score, values_fill = 0)

parent_type <- function(x) {
  if (grepl("^Epi", x)) "Epithelial" else if (grepl("_EC$|iEC$|cycEC$|circEC$", x)) "Endothelial" else if (grepl("CAF|Fibroblast", x)) "Fibroblast" else if (x %in% c("Pericyte", "SMC")) "Mural" else if (grepl("Macro|Mono|cDC|pDC|Mast", x)) "Myeloid" else if (grepl("_B$|^B$|Plasma|Plasmablast|Bgc|Bmem|^Bn$|^ABC$|NR4A2", x)) "Bcell" else "Tcell"
}
score$global_celltype <- vapply(score$name, parent_type, character(1))
score <- score %>% mutate(boundary_advantage = Bdy - (Mal + Normal) / 2) %>% arrange(desc(boundary_advantage))
write.csv(score, file.path(out_dir, "Supplementary_Fig18_betweenness_centrality.csv"), row.names = FALSE)
write.csv(filter(score, Bdy > Mal, Bdy > Normal), file.path(out_dir, "Supplementary_Fig18_boundary_enriched_subtypes.csv"), row.names = FALSE)

# "Reproduce plot" ------------------------------------------------------------
long <- score %>% pivot_longer(c(Mal, Bdy, Normal), names_to = "Category", values_to = "Score")
long$Category <- factor(long$Category, levels = c("Mal", "Bdy", "Normal"))
colors <- c(Mal = "#FF0000", Bdy = "#008B00", Normal = "#4169E1")
types <- c("Epithelial", "Tcell", "Bcell", "Myeloid", "Endothelial", "Fibroblast", "Mural")
plots <- lapply(types, function(type) {
  dat <- filter(long, global_celltype == type)
  ordering <- score %>% filter(global_celltype == type) %>% arrange(boundary_advantage) %>% pull(name)
  ggplot(dat, aes(factor(name, levels = ordering), Score, fill = Category)) +
    geom_col(position = "dodge", color = "black", linewidth = 0.15) + coord_flip() +
    scale_fill_manual(values = colors) + theme_bw() +
    labs(title = type, x = NULL, y = "Betweenness Centrality") +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), axis.text = element_text(size = 6), axis.title = element_text(size = 8), legend.position = "none")
})
pdf(file.path(plot_dir, "Supplementary_Fig18.pdf"), width = 8, height = 7)
for (plot in plots) print(plot)
dev.off()

message("Supplementary Figure 18 complete: ", dirname(out_dir))
