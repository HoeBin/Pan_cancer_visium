"""Extended_Fig2 - Ro/e heatmaps of cell subtypes across cancer types (panel e).

For every major cell class (global_anno_level1) the observed cell count of
each sub_anno x Cancer_type pair is compared with the count expected under
independence (Pearson chi-square contingency table). Ro/e = observed /
expected: > 1 enrichment, < 1 depletion, annotated with the symbols of
Zheng et al. 2021 (+++ > 3, ++ > 1, + >= 0.2, ± > 0, − = 0). One heatmap
per major cell class on a coolwarm diverging scale centred at Ro/e = 1
(white); every cell shows the Ro/e value and its symbol.

Input  : the single-cell h5ad (H5AD below; obs columns global_anno_level1,
         sub_anno, Cancer_type, Data, Sample - X is never read)
Outputs: ../figures/roe_<global_anno_level1>_sub_anno_by_Cancer_type.pdf
         ../figures/roe_sub_anno_by_Cancer_type.xlsx
             (Ro/e matrix, one sheet per major cell class)
         ../figures/celltype_cancertype_ncells_nsamples_ndata.xlsx
             (cell / unique-sample / unique-dataset counts, one sheet
              per major cell class, three blocks per sheet)
"""
import warnings
warnings.filterwarnings('ignore')
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import chi2_contingency

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42     # keep text editable in the PDF
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'Arial:italic'
plt.rcParams['mathtext.rm'] = 'Arial'

FS_TITLE = 22       # figure title (global_anno_level1)
FS_XLABEL = 15      # cancer-type tick labels
FS_YLABEL = 14      # sub_anno tick labels
FS_CELL = 12        # Ro/e value + symbol inside each cell
FS_CBLABEL = 16     # colorbar label
FS_CBTICK = 13      # colorbar tick numbers
VMAX_MIN = 2.0      # colour scale top = max(VMAX_MIN, 99th pct of Ro/e)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                      # ../ = Sup_Fig6
FIGDIR = os.path.join(ROOT, 'figures')
H5AD = ('/path/to/h5ad_annotated/submit_figure_input/'
        'h5ad/_sc_ALL_ligh_260730.h5ad')

# major cell classes in panel order
GROUP_ORDER = ['Epithelial', 'T cell', 'B cell', 'Myeloid',
               'Endothelial', 'Fibroblast', 'Mural']
OBS_COLS = ['global_anno_level1', 'sub_anno', 'Cancer_type', 'Data',
            'Sample']


def roe_symbol(v):
    """Enrichment symbol shown in the heatmap cells (Zheng et al. 2021)."""
    if v > 3:
        return '+++'
    if v > 1:
        return '++'
    if v >= 0.2:
        return '+'
    if v > 0:
        return '±'
    return '−'


def load_obs(path):
    """obs of the h5ad as str columns (backed mode: X is not loaded)."""
    adata = ad.read_h5ad(path, backed='r')
    try:
        obs = adata.obs[OBS_COLS].astype(str).copy()
    finally:
        adata.file.close()
    return obs


def count_tables(obs, cancers):
    """sub_anno x Cancer_type tables: cells, unique samples, unique datasets."""
    n_cells = pd.crosstab(obs['sub_anno'], obs['Cancer_type'])
    n_cells = n_cells[[c for c in cancers if c in n_cells.columns]]

    def n_unique(col):
        return (obs.groupby(['sub_anno', 'Cancer_type'])[col].nunique()
                .unstack(fill_value=0)
                .reindex(index=n_cells.index, columns=n_cells.columns,
                         fill_value=0))

    return n_cells, n_unique('Sample'), n_unique('Data')


def roe_table(n_cells):
    """Ro/e = observed / chi-square expected count; also the test result."""
    chi2, p, dof, expected = chi2_contingency(n_cells)
    roe = pd.DataFrame(n_cells.values / expected,
                       index=n_cells.index, columns=n_cells.columns)
    return roe, chi2, p, dof


def plot_heatmap(gname, roe, out_pdf, png_copy=None):
    """roe: sub_anno x Cancer_type DataFrame."""
    n_r, n_c = roe.shape
    fig, ax = plt.subplots(figsize=(0.7 * n_c + 3, 0.45 * n_r + 2))
    vmax = max(VMAX_MIN, np.percentile(roe.values, 99))
    norm = TwoSlopeNorm(vmin=0, vcenter=1, vmax=vmax)   # white = Ro/e 1
    im = ax.imshow(roe.values, cmap='coolwarm', norm=norm, aspect='auto')
    for i in range(n_r):
        for j in range(n_c):
            v = roe.iloc[i, j]
            # white text on strongly coloured cells
            tc = 'white' if (norm(v) > 0.85 or norm(v) < 0.15) else 'black'
            ax.text(j, i, f'{v:.2f}\n{roe_symbol(v)}', ha='center',
                    va='center', fontsize=FS_CELL, color=tc)
    ax.set_xticks(np.arange(n_c))
    ax.set_xticklabels(roe.columns, fontsize=FS_XLABEL, rotation=45,
                       ha='right', rotation_mode='anchor')
    ax.set_yticks(np.arange(n_r))
    ax.set_yticklabels(roe.index, fontsize=FS_YLABEL)
    ax.set_title(gname, fontsize=FS_TITLE)
    ax.grid(False)
    cax = ax.inset_axes([1.005, 0.0, 0.04, 0.4])   # [x, y, w, h] axes coords
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=FS_CBTICK)
    cb.set_label('Ro/e', fontsize=FS_CBLABEL)
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches='tight')
    if png_copy:
        fig.savefig(png_copy, dpi=110, bbox_inches='tight')
    plt.close(fig)
    print(f'saved {out_pdf}')


if __name__ == '__main__':
    os.makedirs(FIGDIR, exist_ok=True)
    png_dir = os.environ.get('VIOLIN_PNG_DIR')
    obs_all = load_obs(H5AD)
    print('cells per global_anno_level1:')
    print(obs_all['global_anno_level1'].value_counts().to_string())
    cancers = sorted(obs_all['Cancer_type'].unique())   # shared x-axis order

    roe_dict, count_dict = {}, {}
    for gname in GROUP_ORDER:
        obs = obs_all[obs_all['global_anno_level1'] == gname]
        if obs.empty:
            continue
        n_cells, n_samples, n_data = count_tables(obs, cancers)
        roe, chi2, p, dof = roe_table(n_cells)
        roe_dict[gname] = roe
        count_dict[gname] = (n_cells, n_samples, n_data)
        print(f'{gname}: {n_cells.shape[0]} sub_annos x {n_cells.shape[1]} '
              f'cancer types, chi-square p = {p:.2e} (df = {dof})')
        safe = gname.replace(' ', '_')                  # T cell -> T_cell
        plot_heatmap(gname, roe,
                     os.path.join(FIGDIR,
                                  f'roe_{safe}_sub_anno_by_Cancer_type.pdf'),
                     png_copy=png_dir and os.path.join(
                         png_dir, f'sf6_{safe}_roe.png'))

    # Ro/e matrices, one sheet per major cell class
    xlsx = os.path.join(FIGDIR, 'roe_sub_anno_by_Cancer_type.xlsx')
    with pd.ExcelWriter(xlsx, engine='openpyxl') as writer:
        for gname, roe in roe_dict.items():
            roe.round(3).to_excel(writer, sheet_name=gname)
    print(f'saved {xlsx}')

    # cell / unique-sample / unique-dataset counts, three blocks per sheet
    xlsx = os.path.join(FIGDIR,
                        'celltype_cancertype_ncells_nsamples_ndata.xlsx')
    with pd.ExcelWriter(xlsx, engine='openpyxl') as writer:
        for gname, (n_cells, n_samples, n_data) in count_dict.items():
            n_cells.to_excel(writer, sheet_name=gname, startrow=1)
            ws = writer.sheets[gname]
            ws.cell(row=1, column=1, value='Cell count')
            start2 = len(n_cells) + 4
            n_samples.to_excel(writer, sheet_name=gname, startrow=start2)
            ws.cell(row=start2, column=1, value='N samples (unique Sample)')
            start3 = start2 + len(n_samples) + 3
            n_data.to_excel(writer, sheet_name=gname, startrow=start3)
            ws.cell(row=start3, column=1, value='N datasets (unique Data)')
    print(f'saved {xlsx}')
