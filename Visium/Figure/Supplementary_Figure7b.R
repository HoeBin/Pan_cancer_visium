# 분석 목적
#   - scRNA-seq H5AD의 저장된 UMAP 좌표에서 Supplementary Fig.7b를 재현한다.
#
# 분석 흐름
#   1. H5AD에서 expression matrix를 제외하고 UMAP 및 암종 annotation만 읽는다.
#   2. 전체 UMAP을 회색 배경으로 표시하고 암종별 2D density contour를 계산한다.
#   3. 암종별 세포 수와 UMAP 중심·분산을 CSV로 저장한다.
#
# 주요 출력
#   - Supplementary_Fig7b.pdf
#   - 암종별 UMAP 요약 CSV
#
# 출력 위치
#   - Figure/Output/Supplementary_Figure7/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({library(hdf5r); library(data.table); library(dplyr); library(ggplot2)})
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg)); root_dir <- dirname(script_path)
source(file.path(root_dir, "common/config.R"))
out_dir <- file.path(root_dir, "Output/Supplementary_Figure7/Tables"); plot_dir <- file.path(root_dir, "Output/Supplementary_Figure7/Plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
h5ad_file <- file.path(PROCESSED_DATA_DIR, "h5ad_all/_sc_ALL_260730.h5ad")

read_categorical <- function(handle, path) {
  categories <- handle[[paste0(path, "/categories")]][]
  codes <- as.integer(handle[[paste0(path, "/codes")]][])
  out <- rep(NA_character_, length(codes)); keep <- codes >= 0
  out[keep] <- categories[codes[keep] + 1L]; out
}

# "Read reduced dimensions only" ---------------------------------------------
h5 <- H5File$new(h5ad_file, mode = "r")
umap <- h5[["obsm/X_umap"]]$read()
organ <- read_categorical(h5, "obs/Organ_orig")
h5$close_all()
if (nrow(umap) == 2) umap <- t(umap)
stopifnot(nrow(umap) == length(organ))

organ_map <- c(Breast = "BRCA", Colorectal = "COCA", Endometrium = "UECA", `Head and Neck` = "HNCA", Kidney = "KICA", Liver = "LICA", Lung = "LUCA", Ovary = "OVCA", Pancreas = "PACA", Prostate = "PRCA", Skin = "SKCA", Stomach = "STCA", Thyroid = "THCA")
dat <- data.table(UMAP1 = umap[, 1], UMAP2 = umap[, 2], cancer_type = unname(organ_map[organ]))
dat <- dat[!is.na(cancer_type)]
summary_table <- dat[, .(n_cells = .N, mean_UMAP1 = mean(UMAP1), mean_UMAP2 = mean(UMAP2), sd_UMAP1 = sd(UMAP1), sd_UMAP2 = sd(UMAP2)), by = cancer_type][order(cancer_type)]
write.csv(summary_table, file.path(out_dir, "Supplementary_Fig7b_key_results_umap_by_cancer.csv"), row.names = FALSE)

# Limit background points deterministically for PDF size; density uses all cells.
set.seed(1)
background <- dat[sample.int(.N, min(.N, 250000))]
cancer_order <- c("BRCA", "COCA", "UECA", "HNCA", "KICA", "LICA", "LUCA", "OVCA", "PACA", "PRCA", "SKCA", "STCA", "THCA")
dat$cancer_type <- factor(dat$cancer_type, levels = cancer_order)
p <- ggplot() +
  geom_point(data = background, aes(UMAP1, UMAP2), color = "grey85", size = 0.08, alpha = 0.5) +
  stat_density_2d(data = dat, aes(UMAP1, UMAP2, fill = after_stat(level)), geom = "polygon", alpha = 0.75, bins = 5, color = "black", linewidth = 0.25) +
  scale_fill_gradient(low = "#FEE5D9", high = "#CB181D") + facet_wrap(~cancer_type, ncol = 5) +
  coord_equal() + theme_void() + theme(strip.text = element_text(size = 12, face = "bold"), legend.position = "none")
pdf(file.path(plot_dir, "Supplementary_Fig7b.pdf"), width = 14, height = 8, useDingbats = FALSE)
print(p)
dev.off()
message("Supplementary Figure 7b complete: ", dirname(out_dir))
