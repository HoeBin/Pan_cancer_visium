# 분석 목적
#   - Source_Data_Fig2.xlsx로부터 Figure 2e, 2g, 2h를 재현한다.
#
# 분석 흐름
#   1. Fig.2e left/right(scRNA-seq, Visium cell-type proportion)를 암종별 boxplot으로 그린다.
#   2. Fig.2g(compartment별 marker gene DE 결과)에서 상위 유전자를 dot plot으로 그린다.
#   3. Fig.2h(Global/Myeloid compartment proportion)를 heatmap으로 그린다.
#
# 주요 출력
#   - Fig2e_scRNAseq.pdf, Fig2e_Visium.pdf
#   - Fig2g.pdf
#   - Fig2h.pdf
#
# 출력 위치
#   - Output/Figure2/
#
# 참고: Fig.2a-d(scRNA-seq/Visium UMAP), Fig.2f(malignancy UMAP density)는 Source Data가
# 제공되지 않는 도식/embedding plot이라 재현 대상에서 제외한다. Fig.2g는 Source Data에
# compartment 자기 자신에 대한 DE 통계(scores/logFC/pvals)만 있고, scanpy dotplot 특유의
# 3-group 교차 발현값(다른 compartment에서의 fraction/mean expression)은 없으므로,
# dot 크기=|logFC|, 색=score로 표현하는 근사 dot plot으로 재현한다.

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(ComplexHeatmap); library(circlize); library(viridisLite); library(grid)
})
script_arg  <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(script_path), "common/config.R"))

out_dir  <- file.path(output_dir, "Figure2/Tables")
plot_dir <- file.path(output_dir, "Figure2/Plots")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# "Figure 2e" -------------------------------------------------------------
sc  <- read_source("Source_Data_Fig2.xlsx", "Fig.2e left")
vis <- read_source("Source_Data_Fig2.xlsx", "Fig.2e right")
names(sc)  <- c("cancer_type", "sample_id", "celltype", "count", "proportion")
names(vis) <- c("cancer_type", "celltype", "sample_id", "proportion")
sc$celltype  <- recode(sc$celltype,  `T cell` = "Tcell", `B cell` = "Bcell")
vis$celltype <- recode(vis$celltype, `T cell` = "Tcell", `B cell` = "Bcell")
sc  <- sc[sc$cancer_type %in% cancer_order, ]     # nonESCA 논문 버전과 동일하게 ESCA 제외
vis <- vis[vis$cancer_type %in% cancer_order, ]
sc$celltype   <- factor(sc$celltype,  levels = celltype_order)
vis$celltype  <- factor(vis$celltype, levels = celltype_order)
sc$cancer_type  <- factor(sc$cancer_type,  levels = cancer_order[cancer_order %in% sc$cancer_type])
vis$cancer_type <- factor(vis$cancer_type, levels = cancer_order[cancer_order %in% vis$cancer_type])
write.csv(sc,  file.path(out_dir, "Fig2e_scRNAseq_plot_data.csv"), row.names = FALSE)
write.csv(vis, file.path(out_dir, "Fig2e_Visium_plot_data.csv"),   row.names = FALSE)

plot_prop <- function(dat, title) {
  ggplot(dat, aes(celltype, proportion)) +
    geom_jitter(size = 0.5, width = 0.2) +
    geom_boxplot(aes(fill = celltype), outlier.shape = NA, alpha = 0.8) +
    facet_wrap(~cancer_type, ncol = 5, scales = "free_y", drop = TRUE) +
    scale_fill_manual(values = celltype_colors, drop = FALSE) +
    scale_x_discrete(labels = celltype_labels, drop = FALSE) +
    base_theme +
    labs(title = title, x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 10))
}
pdf(file.path(plot_dir, "Fig2e_scRNAseq.pdf"), width = 9, height = 6)
print(plot_prop(sc, "scRNA-seq"))
dev.off()
pdf(file.path(plot_dir, "Fig2e_Visium.pdf"), width = 9, height = 6)
print(plot_prop(vis, "Visium"))
dev.off()

# "Figure 2g" -------------------------------------------------------------
deg <- read_source("Source_Data_Fig2.xlsx", "Fig.2g")
deg <- deg[!is.na(deg$Gene), ]
deg$Compartment <- factor(deg$Compartment, levels = c("Normal", "Boundary", "Malignant"))
top_n <- 10
top_genes <- deg %>%
  group_by(Compartment) %>%
  arrange(desc(scores), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()
gene_order <- top_genes %>% arrange(Compartment) %>% pull(Gene) %>% unique()
top_genes$Gene <- factor(top_genes$Gene, levels = gene_order)
write.csv(top_genes, file.path(out_dir, "Fig2g_top_marker_genes.csv"), row.names = FALSE)

p2g <- ggplot(top_genes, aes(x = Gene, y = Compartment)) +
  geom_point(aes(size = abs(logFC), color = scores)) +
  scale_color_gradient(low = "#FDE0DD", high = "#8B0000", name = "score") +
  scale_size_continuous(range = c(2, 8), name = "|logFC|") +
  base_theme +
  labs(title = "Figure 2g: top compartment marker genes (own-group DE stats)", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11))
pdf(file.path(plot_dir, "Fig2g.pdf"), width = 12, height = 4)
print(p2g)
dev.off()

# "Figure 2h" -------------------------------------------------------------
h <- read_source("Source_Data_Fig2.xlsx", "Fig.2h")
h <- h[!is.na(h$CellType), ]
col_fun <- colorRamp2(c(0, 0.9, 1), viridis(3, option = "C"))
heatmaps <- list()
for (ctype in unique(h$Type)) {
  x <- h %>% filter(Type == ctype) %>% select(CellType, Compartment, ScaledMedian) %>%
    pivot_wider(names_from = Compartment, values_from = ScaledMedian) %>% as.data.frame()
  rownames(x) <- x$CellType; x$CellType <- NULL
  mat <- as.matrix(x[, intersect(mal_order, colnames(x)), drop = FALSE])
  heatmaps[[ctype]] <- Heatmap(mat, name = "proportion", column_title = ctype,
    col = col_fun, cluster_columns = FALSE, cluster_rows = TRUE,
    column_names_rot = 45, row_names_gp = gpar(fontsize = 15),
    column_names_gp = gpar(fontsize = 15),
    column_title_gp = gpar(fontsize = 20, fontface = "bold"),
    heatmap_legend_param = list(at = c(0, 1), labels = c("min", "max"),
      title_gp = gpar(fontsize = 15, fontface = "bold"),
      labels_gp = gpar(fontsize = 12)))
}
pdf(file.path(plot_dir, "Fig2h.pdf"), width = 7, height = 7)
for (ctype in names(heatmaps)) draw(heatmaps[[ctype]], padding = unit(c(2, 40, 2, 2), "mm"))
dev.off()

message("Figure 2 complete: ", dirname(out_dir))
