# Figure reproduction (Main Figures 1, 2, 4, 5 + Extended Data Figures 4, 5, 6)

이 폴더는 논문에 공식 배포된 `Source_Data/*.xlsx`만을 입력으로 사용해 Main/Extended
Data Figure의 데이터 기반 패널을 재현하는 self-contained R 코드다. 저장소를 clone한
뒤 별도 경로 설정 없이 바로 실행할 수 있다.

이전 버전(원시 분석 산출물을 `common/config.R`의 `PROJECT_DIR`/`DECON_DIR` 등으로
가리켜 재현하던 스크립트: `Figure2e_2h.R`, `Figure3.R`, `Figure4.R`, `Figure5a_5d.R`,
`Supplementary_Figure*.R`)를 이 Source-Data 기반 스크립트로 교체했다. 그 결과 Figure3
전체, Figure4d/f/g, Figure5d, Supplementary Figure 패널의 재현 코드는 이 폴더에서
빠졌다 — 필요하면 git 히스토리에서 복구할 수 있다.

## 디렉토리 구조

```
Figure/
  common/config.R
  Figure1.R, Figure2.R, Figure4.R, Figure5.R
  ExtendedFigure4.R, ExtendedFigure5.R, ExtendedFigure6.R
  Source_Data/*.xlsx                     — 논문 공식 배포 Source Data (입력)
  Output/FigureN/{Plots,Tables}/         — Main Figure 재현 PDF 및 CSV
  Output/ExtendedFigureN/{Plots,Tables}/ — Extended Data Figure 재현 PDF 및 CSV
```

Extended Fig.1, 2는 로컬에는 `ExtendedFigure1.R`/`ExtendedFigure2.R`로 존재하지만
`.gitignore`에 등록되어 이 GitHub 저장소에는 올라가지 않는다.

## 실행

```bash
cd Visium/Figure
Rscript Figure1.R                   # 각 Figure는 독립적으로 실행 가능
Rscript ExtendedFigure4.R
```

## 범위 및 한계

### Main Figure

Source Data가 제공된 패널만 재현했다. 아래 패널은 개별 spot 좌표, UMAP embedding,
BioRender 일러스트 등 원본 이미지/도식 데이터가 Source Data에 없어 제외했다:

- Figure 1b (BioRender workflow schematic)
- Figure 2a-d, 2f (scRNA-seq/Visium UMAP embedding, marker dot plot, malignancy UMAP density)
- Figure 4e (PACA 샘플 spatial multi-panel)
- Figure 5c (UECA/HNCA/BRCA/OVCA/LUCA 샘플 spatial co-enrichment 이미지)

아래 패널은 Source Data가 있음에도 이번 범위에서 제외했다:

- Figure 3 전체 (a-h)
- Figure 4d (HSPA6+ Macro-tCAF LR pair rank), 4f (hazard ratio forest plot), 4g (Kaplan-Meier)
- Figure 5d (scRNA-seq/Visium LR strength 암종별 density plot 및 상관관계)

일부 패널은 Source Data의 정보량 제약으로 완전히 동일하지 않고 근사적으로만
재현했다:

- **Fig2g**: Source Data는 각 유전자가 자신의 compartment 그룹 내에서 얻은 DE 통계
  (scores/logFC/pvals)만 제공하고, scanpy dotplot 특유의 3-group 교차 발현값(다른
  compartment에서의 발현 비율·평균)은 없다. compartment별 상위 10개 유전자를
  dot 크기=|logFC|, 색=score로 표현하는 근사 dot plot으로 재현했다.

나머지 패널(Fig1a, Fig2e/h, Fig4a/b/c, Fig5a/b)은 Source Data의 값과 통계량을 그대로
사용해 재현했다.

### Extended Data Figure

Extended Fig.1, 2, 3, 7, 8, 10은 이번 GitHub 범위에서 제외했다.

- **Extended Fig.1, 2**: 재현 코드는 로컬에 있지만(위 참고) 이번 저장소 공개 범위에서
  제외했다.
- **Extended Fig.3, 7, 8, 10**: Source Data는 존재하나 이번 범위에서 제외했다.
- **Extended Fig.4**: a-d(compartment별 CNV score/spot 수 boxplot, global/subtype
  proportion과 CNV score의 상관관계 bar plot) 전체를 재현했다.
- **Extended Fig.5**: a(major cell-type pair의 compartment별 비교, UpSet-matrix
  스타일 pair indicator 포함), b(Epi-Epi/Epi-TME/TME-TME compartment 비교) 전체를
  재현했다.
- **Extended Fig.6**: Malignant/Normal co-enrichment network를 재현했다. 원본은
  Cytoscape 수작업 배치라 노드 좌표는 다르지만(Source Data의 edge list로
  force-directed layout을 새로 계산), 네트워크 구조(연결 관계)는 동일하다. subtype
  이름에 그리스 문자(gamma-delta T)가 있어 `cairo_pdf()` 장치를 사용했다(기본
  `pdf()`는 해당 글자를 렌더링하지 못함).

Extended Fig.5a의 pair 순서와 Extended Fig4d·Fig6의 subtype -> lineage(색상) 매핑은
`common/config.R`의 `celltype_order`, `load_subtype_lineage()`를 공유해서 쓴다.

## 다음 단계

- Figure 3, Fig4d/f/g, Fig5d / Extended Fig.3, 7, 8, 10은 Source Data가 있으므로
  필요 시 언제든 같은 패턴으로 다시 추가할 수 있다.
- Supplementary Figure 패널(공식 Source Data에 데이터 존재)도 동일한 패턴
  (`common/config.R` 재사용, `SupplementaryFigureN.R` 스크립트 추가)으로 이어서
  작성할 수 있다.
