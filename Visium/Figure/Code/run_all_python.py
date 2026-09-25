#!/usr/bin/env python3
"""Paper_info/Source_Data의 공식 Source Data만을 입력으로 사용해 Figure3 전체, Fig.4(a,b,d,f,g),
Fig.5d, Ext.Fig.3/5b/6/7/8/10, Supp.Fig.11-13의 Python 패널을 순서대로 실행한다. Fig.4c, Fig.5(a,b)와
Ext.Fig.1/2/4는 R로 재현되며 run_all.R / run_all_extended.R이 담당한다 (겹치는 패널은 이 스크립트,
즉 Fig345_code_sourcedata_upload_최종_use_Ver2에서 가져온 Python 코드가 우선).

개별 스크립트도 `python Code/Figure3/*.py` 등으로 독립 실행 가능하며 서로 의존하지 않는다
(Fig4e_ExtDataFig9_slide_panels.py / ExtDataFig9_grid_3x3.py 제외: 후자가 전자의 결과물을 사용).

Fig.4e/Ext.Fig.9와 Fig.5c는 공식 배포된 Source Data가 아니라 원본 Visium h5ad가 있어야 그려진다.
PCASSO_VISIUM_H5AD 환경변수가 설정되어 있을 때만 실행하고, 없으면 건너뛴다.
"""
import os
import subprocess
import sys
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent

SOURCE_DATA_SCRIPTS = [
    "Figure3/Fig3d_ExtDataFig6_8_networks.py",
    "Figure3/Fig3ef_ExtDataFig7ef_centrality.py",
    "Figure3/Fig3cgh_ExtDataFig7cgh.py",
    "Figure3/Fig3b_ExtDataFig7b.py",
    "Figure3/ExtDataFig5b_7d_pair_categories.py",
    "Figure3/Suppfig11_bc_by_compartment.py",
    "Figure4/Fig4ab_signature_panels.py",
    "Figure4/Fig4d_ranked_scatter.py",
    "Figure4/Fig4fg_SupFig12-13_survival.py",
    "Figure4/ExtDataFig10_bulk_correlation.py",
    "Figure5/Fig5d_sc_vs_visium_scatter.py",
    "ExtendedFigure3.py",
]

H5AD_SCRIPTS = [
    "Figure4/Fig4e_ExtDataFig9_slide_panels.py",
    "Figure4/ExtDataFig9_grid_3x3.py",
    "Figure5/Fig5c_spatial_lr_overlap.py",
]


def run(rel_path: str) -> None:
    print(f"==== Running {rel_path} ====")
    subprocess.run([sys.executable, str(CODE_DIR / rel_path)], check=True, cwd=CODE_DIR)


if __name__ == "__main__":
    for rel_path in SOURCE_DATA_SCRIPTS:
        run(rel_path)

    if os.environ.get("PCASSO_VISIUM_H5AD"):
        for rel_path in H5AD_SCRIPTS:
            run(rel_path)
    else:
        print("PCASSO_VISIUM_H5AD not set: skipping Fig.4e/ExtFig.9/Fig.5c "
              f"({', '.join(H5AD_SCRIPTS)}) — these need the deposited Visium h5ad, not Source Data.")
