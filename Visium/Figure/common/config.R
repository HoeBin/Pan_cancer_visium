# Figure 재현에 필요한 처리된 데이터 경로 설정.
# 아래 각 경로를 실제 데이터 위치로 수정한 뒤 Figure 스크립트를 실행하세요.

# 분석 파이프라인 산출물 루트: Cell2Location_SubNum4/ 이하에 celltype mapping,
# co-enrichment network, malignancy/CNV 분석, Visium gene-set score,
# Fig5 관련 supplementary table 등을 포함
PROJECT_DIR <- "/path/to/input/cell2loc_project"

# cell2location deconvolution 결과 루트: 암종별 cell2location_map_*_30000epoch/
# 서브디렉토리에 Obs.csv, Obs_abundance.csv 포함
DECON_DIR <- "/path/to/input/decon_results"

# 전처리된 부가 데이터: malignancy/CNV 영역 주석(region_output_2025_02_02.csv),
# Visium gene-set score(Visium_Scoring_251214.csv), scRNA-seq h5ad(h5ad_all/)
PROCESSED_DATA_DIR <- "/path/to/input/processed_data"

# STopover 공동국소화(co-localization) 결과: global/, subtype/ 서브디렉토리
STOPOVER_DIR <- "/path/to/input/stopover_results"

# Ligand-receptor 공동국소화(Jaccard) 결과
LR_DIR <- "/path/to/input/lr_results"
