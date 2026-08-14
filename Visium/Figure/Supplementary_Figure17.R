# 분석 목적
#   - Supplementary Figure 17(Malignant·Normal 세포 유형 co-enrichment network)의 입력 table을 정리한다.
#
# 분석 흐름
#   1. celltype_mapping.csv로 subtype label을 최종 명칭으로 변환할 rename key를 만든다.
#   2. Malignant·Normal 영역별 network edge(Jaccard co-localization score > 0.5 threshold)/node table을 불러온다.
#   3. label을 정리해 CSV로 저장한다.
#
# 주요 출력
#   - Supplementary_Fig17a(Malignant)/17b(Normal) network edge/node CSV
#   - Cytoscape에서 수기로 배치한 다이어그램이라 PDF 재현 대상은 아님
#
# 출력 위치
#   - Figure/Output/Supplementary_Figure17/

# "Libraries and paths" -------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
root_dir <- dirname(script_path)
source(file.path(root_dir, "common/config.R"))
project <- PROJECT_DIR

out_dir <- file.path(root_dir, "Output/Supplementary_Figure17/Tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# "Celltype label mapping" -----------------------------------------------------
celltype_map <- fread(file.path(project, "Cell2Location_SubNum4/Data_Info/celltype_mapping.csv"))
rename_key <- setNames(celltype_map$after_celltype, celltype_map$before_celltype)

# "Panels a-b: Malignant and Normal co-enrichment networks" -------------------
node_dir <- file.path(project, "Cell2Location_SubNum4/Network/Stopover2")
regions <- c(Mal = "a_malignant", Normal = "b_normal")
for (region in names(regions)) {
  suffix <- regions[[region]]

  edge <- fread(file.path(node_dir, paste0("Cor_cluster_", region, "_cut0.5.csv")))
  edge[, Var1 := fifelse(Var1 %chin% names(rename_key), rename_key[Var1], Var1)]
  edge[, Var2 := fifelse(Var2 %chin% names(rename_key), rename_key[Var2], Var2)]
  edge[, Var3 := paste(Var1, "(interacts with)", Var2)]
  write.csv(edge, file.path(out_dir, paste0("Supplementary_Fig17", suffix, "_network_edges.csv")), row.names = FALSE)

  node <- fread(file.path(node_dir, paste0("Cor_cluster_", region, "_cut0.5.csv_1 default node.csv")))
  node[, name := fifelse(name %chin% names(rename_key), rename_key[name], name)]
  node[, `shared name` := name]
  write.csv(node, file.path(out_dir, paste0("Supplementary_Fig17", suffix, "_network_nodes.csv")), row.names = FALSE)
}

message("Supplementary Figure 17 complete: ", dirname(out_dir))
