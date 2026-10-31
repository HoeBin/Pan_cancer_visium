import os
import scanpy as sc

# 특정 경로 설정

directory_path = '/path/to/h5ad_raw'  # .h5ad 파일들이 있는 디렉토리 경로

# QC를 거친 데이터들을 저장할 리스트 초기화
adatas = []

# 디렉토리 내의 모든 .h5ad 파일에 대해 반복 수행
for filename in os.listdir(directory_path):
    if filename.endswith('.h5ad'):
        # 파일 경로 생성
        file_path = os.path.join(directory_path, filename)
        
        # .h5ad 파일 불러오기
        adata = sc.read_h5ad(file_path)
        adata.obs_names_make_unique()
        
        # 1. mt와 ribo 유전자 설정
        adata.var['mt'] = adata.var_names.str.startswith('MT-') 
        adata.var['ribo'] = adata.var_names.str.startswith(("RPS", "RPL"))
        
        # 2. QC metrics 계산
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        # 3. 총 UMI 수 계산 및 필터링
        sc.pp.filter_cells(adata, min_counts=1000)  # 1000 UMI 미만 제거
        sc.pp.filter_cells(adata, min_genes=500)    # 500개 미만의 유전자 제거
        
        # 4. 7000개 이상의 유전자가 검출된 세포 제거 (잠재적 이중)
        adata = adata[adata.obs.n_genes_by_counts < 7000, :]
        
        # 5. 미토콘드리아 유전자 비율이 20% 이상인 세포 제거
        adata = adata[adata.obs.pct_counts_mt <= 20, :]

        # QC된 데이터 리스트에 추가
        adatas.append(adata)

# 모든 QC된 데이터들을 하나로 결합 (concat)
combined_adata = sc.concat(adatas, join='inner')

# 최종 결합된 데이터를 저장
output_file = '/path/to/h5ad_raw/concat_data/01_QC_concat.h5ad'
combined_adata.write_h5ad(output_file)

# 결과 확인 (선택 사항)
print(f"최종 결합된 데이터: {combined_adata}")
