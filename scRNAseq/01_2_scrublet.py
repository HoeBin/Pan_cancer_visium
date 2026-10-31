import os
import pandas as pd
import numpy as np
import scanpy as sc
from tqdm import tqdm
import scrublet as scr

def run_scrublet(adata, column):
    # adata에서 20개 이하의 sample을 가진 행을 제거
    sample_counts = adata.obs[column].value_counts()
    valid_samples = sample_counts[sample_counts > 20].index
    adata = adata[adata.obs[column].isin(valid_samples)]
    
    # 결과를 저장할 빈 DataFrame을 생성
    results = pd.DataFrame(index=adata.obs.index)
    
    # 에러가 발생한 샘플을 저장할 리스트
    error_samples = []
    
    unique_samples = adata.obs[column].unique()
    for sample in tqdm(unique_samples, desc='Processing samples'):
        try:
            sample_data = adata[adata.obs[column] == sample].copy()
            
            # Scrublet을 사용하여 doublet scores와 예측값을 계산
            scrub = scr.Scrublet(sample_data.X)
            sample_data.obs['scrublet_doublet_scores'], sample_data.obs['scrublet_predicted_doublets'] = scrub.scrub_doublets(verbose=False)
            sample_data.obs['scrublet_predicted_doublets'] = sample_data.obs['scrublet_predicted_doublets'].astype(str)
            
            # 결과를 저장할 DataFrame에 결과를 저장
            results.loc[sample_data.obs.index, "scrublet_doublet_scores"] = sample_data.obs['scrublet_doublet_scores']
            results.loc[sample_data.obs.index, "scrublet_predicted_doublets"] = sample_data.obs['scrublet_predicted_doublets']
        except Exception as e:
            print(f"Error processing sample {sample}: {e}")
            error_samples.append(sample)

    # 에러가 발생한 샘플 정보 출력
    if error_samples:
        print("Errors occurred with the following samples:")
        for sample in error_samples:
            print(sample)

    return results, error_samples

# 특정 파일 경로 설정
file_path = '/path/to/h5ad_raw/concat_data/01_QC_concat.h5ad'

# 전체 결과를 저장할 데이터프레임 초기화
final_results = pd.DataFrame()

try:
    # .h5ad 파일 불러오기
    adata = sc.read_h5ad(file_path)
    adata.obs_names_make_unique()
    
    # Scrublet 실행
    results, error_samples = run_scrublet(adata, 'Sample')
    
    # 에러가 발생한 경우, 에러 메시지와 파일명을 결과에 추가
    if error_samples:
        error_df = pd.DataFrame({
            'scrublet_doublet_scores': [np.nan] * len(error_samples),
            'scrublet_predicted_doublets': [np.nan] * len(error_samples)
        }, index=error_samples)
        error_df['Error'] = f"Error in file {os.path.basename(file_path)}"
        results = pd.concat([results, error_df])
    
    # 결과를 최종 결과 데이터프레임에 추가
    final_results = pd.concat([final_results, results])

except Exception as e:
    print(f"Error processing file {file_path}: {e}")
    # 에러 발생 시, 해당 파일 이름과 에러 메시지를 기록
    error_df = pd.DataFrame({
        'scrublet_doublet_scores': [np.nan],
        'scrublet_predicted_doublets': [np.nan],
        'Error': [str(e)]
    }, index=[os.path.basename(file_path)])
    final_results = pd.concat([final_results, error_df])

# 최종 결과를 CSV 파일로 저장 (선택 사항)
final_results.to_csv('/path/to/qc_output/scrublet_results_new.csv')

# 최종 결과 출력 (선택 사항)
print(final_results)
