# Figure reproduction package (Main Figures 1-5 + Extended Data Figures 1-10 + Supplementary Fig. 11-13)

이 저장소는 논문에 공식 배포된 `Paper_info/Paper_info/Source_Data/*.xlsx`만을 입력으로
사용해 Main/Extended Data Figure의 데이터 기반 패널을 재현하는 코드(`Code/`)와 그
출력(`Output/`)을 담는다. R과 Python이 섞여 있다 — 각 패널은 원래 검증된 코드베이스의
언어를 그대로 유지했다(아래 "R/Python 혼합" 참고). `Paper_info/Paper_info`는 논문 배포
자료(Source Data, Manuscript, Figure PDF 등) 원본을 그대로 보관하는 용도이며 이 코드가
그 안에 어떤 파일도 쓰지 않는다.

두 코드 계층 중 **Figure layer만** 포함한다: Source Data(xlsx)만 읽어 그림을 그리는
스크립트만 있고, raw h5ad/TCGA counts 등 외부 원본 데이터로부터 통계를 재계산하는
Analysis layer는 포함하지 않는다 — 이 저장소는 공식 Source Data만 입력으로 쓰는
self-contained 재현 패키지다. 예외적으로 h5ad가 있어야만 그려지는 세 이미지 패널
(Fig.4e, Ext.Fig.9, Fig.5c)만 포함하되, 실행 시 외부 Visium h5ad 경로가 필요함을
아래에 표시했다.

## 디렉토리 구조

```
Visium/Figure/
  Code/
    common/config.R, config.py        — 경로, 팔레트, Source Data 리더 (R/Python 짝)
    Data_Info/subtype_lineage_map.csv — Fig.3d/Ext.Fig.6/8 네트워크 노드 색상용 (93 subtype -> 7 lineage)
    Figure1.R, Figure2.R                       — R, 단독
    Figure3/Fig3a_ExtDataFig7a_correlation_heatmap.R — R, Fig.3a + Ext.Fig.7a (Fig.3a는 Input/fig3a_source_data.csv 사용, 아래 "Fig.3a 입력" 참고)
    Figure3/Input/fig3a_source_data.csv        — Fig.3a 상관 행렬 번들 입력(Source Data 아님)
    Figure3/Fig3d_ExtDataFig6_8_networks.py, Fig3ef_ExtDataFig7ef_centrality.py,
           Fig3cgh_ExtDataFig7cgh.py, Fig3b_ExtDataFig7b.py,
           ExtDataFig5b_7d_pair_categories.py, Suppfig11_bc_by_compartment.py
                                                — Python, Fig.3(b-h) + Ext.Fig.5b/6/7/8 + Supp.Fig.11
    Figure4/Fig4ab_signature_panels.py, Fig4d_ranked_scatter.py,
           Fig4fg_SupFig12-13_survival.py, ExtDataFig10_bulk_correlation.py
                                                — Python, Fig.4(a,b,d), Fig.4(f,g)+Supp.Fig.12-13, Ext.Fig.10
    Figure4/Fig4c_pEMT_abundance_correlation.R — R, Fig.4c
    Figure4/Fig4e_ExtDataFig9_slide_panels.py,
           ExtDataFig9_grid_3x3.py  ⚠ 외부 h5ad 필요, 후자는 전자의 결과물을 사용
                                                — Fig.4e, Ext.Fig.9
    Figure4/Input/pemt_scores.csv.gz            — 위 h5ad 스크립트가 쓰는 번들 입력(Source Data 아님)
    Figure5/Fig5ab_pancancer_LR_zscore.R        — R, Fig.5a-b
    Figure5/Fig5d_sc_vs_visium_scatter.py       — Python, Fig.5d
    Figure5/Fig5c_spatial_lr_overlap.py  ⚠ 외부 h5ad 필요 — Fig.5c
    Figure5/Input/selected_slides.csv           — 위 h5ad 스크립트가 쓰는 번들 입력(Source Data 아님)
    ExtendedFigure1.R, ExtendedFigure2.R, ExtendedFigure4.R, ExtendedFigure5.R(5a만) — R
    ExtendedFigure3.py                          — Python, Fig.2의 companion(Ext.Fig.3a-h)
    run_all.R            — R 패널 순서 실행 (Fig1, 2, 4c, 5a-b)
    run_all_extended.R   — R Extended 패널 순서 실행 (ExtFig 1, 2, 4, 5a)
    run_all_python.py    — Python 패널 순서 실행 (Fig.3 전체, Fig.4(a,b,d,f,g), Fig.5d, Ext.Fig.3/5b/6/7/8/10,
                            Supp.11-13; PCASSO_VISIUM_H5AD 설정 시 h5ad 패널도 이어서 실행)
  Output/FigureN/{Plots,Tables}/, ExtendedFigureN/{Plots,Tables}/, SupplementaryFigureN/{Plots,Tables}/
  Paper_info/Paper_info/         — 논문 공식 배포 자료 (입력, read-only)
    Source_Data/*.xlsx           — 이 저장소에는 Source_Data만 포함한다 (원고·Figure PDF/PPTX 등은 제외)
    Figure/, Manuscript.docx, ...
```

## R/Python 혼합, 그리고 겹치는 패널의 우선순위

Fig.3/4/5는 두 개의 원본 코드베이스(R 기반 `Github_code`, Python 기반
`Fig345_code_sourcedata_upload_최종_use_Ver2`)가 부분적으로 겹쳐서 만들어졌다. 겹치는
패널은 2026-09-07 자로 원본 파이프라인과 픽셀 단위까지 검증된 Python 코드를
우선했다:

| 패널 | 우선 채택 | 비고 |
|---|---|---|
| Fig.4a, 4b | Python (`Figure4/Fig4ab_signature_panels.py`) | 같은 268 슬라이드가 3개 compartment에 매칭되는 paired 데이터라 paired Wilcoxon(+Friedman)을 쓴다. 이전 R 버전은 `stat_compare_means`에 `paired=TRUE`가 없어 unpaired Wilcoxon으로 계산되고 있었다. |
| Ext.Fig.5b | Python (`Figure3/ExtDataFig5b_7d_pair_categories.py`) | 같은 이유로 paired Wilcoxon. R 버전(Ext.Fig.5b)은 같은 방식으로 unpaired였다. |
| Ext.Fig.6 (Malignant/Normal co-enrichment network) | Python (`Figure3/Fig3d_ExtDataFig6_8_networks.py`) | Fig.3d/Ext.Fig.8과 같은 스크립트에서 excel 15-digit 정밀도까지 맞춰 노드 배치를 원본과 동일하게 재현한다. |

겹치지 않아 원래 코드를 그대로 쓴 것: Fig.4c(R, `Figure4/Fig4c_pEMT_abundance_correlation.R`), Fig.5a-b(R,
`Figure5/Fig5ab_pancancer_LR_zscore.R`), Ext.Fig.5a(R, `ExtendedFigure5.R`), Fig.1/2 및 Ext.Fig.1/2/4(R, 전부
Dir1이 다루지 않음).

## 범위 및 한계 (Source Data가 없어 제외한 패널)

- Figure 1b (BioRender workflow schematic)
- Figure 2a-d, 2f (scRNA-seq/Visium UMAP embedding, marker dot plot, malignancy UMAP density)
- Fig2g는 Source Data의 정보량 제약으로 근사 dot plot으로만 재현했다(compartment별 상위 10개 유전자, dot 크기=|logFC|, 색=score).
- Extended Fig.1a-b, 1d-e(scRNA-seq UMAP, marker gene dot plot); 1c/1f(Hallmark pathway heatmap)만 재현.
- Extended Fig.6은 Cytoscape 수작업 배치가 아니라 Source Data의 edge list로 force-directed layout을 다시 계산한 것이라 노드 좌표가 다를 수 있지만 excel 15-digit 정밀도로 원본과 동일한 배치를 재현한다(네트워크 구조는 동일).

## 실행

Python 환경: numpy, pandas, scipy, matplotlib, networkx, seaborn, openpyxl, h5py, lifelines(Fig.4f/g),
adjustText(Ext.Fig.10, Ext.Fig.9), pypdf(Ext.Fig.9 3x3 격자), anndata/scanpy(h5ad 패널만). Arial이 설치되어
있으면 사용하고 없으면 기본 sans-serif로 대체한다.

```bash
cd Visium/Figure
Rscript Code/run_all.R                    # Fig.1, 2, 3a(+Ext.Fig.7a), 4c, 5a-b
Rscript Code/run_all_extended.R           # Ext.Fig.1, 2, 4, 5a
python Code/run_all_python.py             # Fig.3 전체, Fig.4(a,b,d,f,g), Fig.5d, Ext.Fig.3/5b/6/7/8/10, Supp.11-13

# h5ad가 있어야 그려지는 세 이미지 패널(Fig.4e, Ext.Fig.9, Fig.5c) — 선택 사항
export PCASSO_VISIUM_H5AD=/real/path/_visium_ALL_260730.h5ad
python Code/run_all_python.py             # 위 환경변수가 설정되어 있으면 h5ad 패널까지 이어서 실행
```

## Fig.3a 입력

`Source_Data_Fig3.xlsx`의 `Fig.3a` 시트는 옛 proportion 기반 상관 값(-0.62~0.61)이라 원고 서술("Epithelial cells
showed near-zero correlations")과 맞지 않는다. 그래서 Fig.3a는 abundance 기반 상관(ESCA 제외 13개 암종,
0.016~0.833)을 담은 `Code/Figure3/Input/fig3a_source_data.csv`를 기본 입력으로 쓴다. Excel이 갱신되면

```bash
Rscript Code/Figure3/Fig3a_ExtDataFig7a_correlation_heatmap.R --fig3a-source=xlsx
```

로 `Fig.3a` 시트를 바로 사용할 수 있다. 기본(csv) 모드에서는 실행 때마다 CSV와 xlsx 시트의 최대 절대 차이를
출력하므로 갱신 여부를 알 수 있다. Ext.Fig.7a는 `Source_Data_Extended_Fig7.xlsx`의 `ED Fig.7a` 시트(BRCA 30 slide
다운샘플, 이미 abundance 기반)를 그대로 읽는다.

개별 스크립트도 각자 독립 실행 가능하다(`Rscript Code/Figure1.R`, `python Code/Figure3/Fig3d_ExtDataFig6_8_networks.py` 등).
`Figure4/ExtDataFig9_grid_3x3.py`는 `Figure4/Fig4e_ExtDataFig9_slide_panels.py`의 결과물을 사용하므로
그 스크립트를 먼저 실행해야 한다.
