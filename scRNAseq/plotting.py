"""Plotting helpers shared by the notebooks.

`umap_with_arrow_axes` is used by the six 07_* sub-annotation notebooks;
`apply_style` and the three ORA functions by the four 08_*_ORA notebooks.
Each was previously duplicated verbatim in every notebook that used it.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from matplotlib.lines import Line2D


# ---------------------------------------------------------------- UMAP panel (07_*)

def umap_with_arrow_axes(adata, color="sub_anno", ncol=3, savefig=None):
    """UMAP coloured by `color`, arrow axes in place of a frame, legend under the panel.

    Colours are taken from `adata.uns[f"{color}_colors"]`; set that beforehand to control them.
    Returns the legend, so a notebook cell ending in this call shows the same repr as before.
    """
    ax = sc.pl.umap(
        adata,
        color=color,
        frameon=False,
        title="",
        legend_loc=None,   # scanpy 기본 레전드 끄기
        show=False
    )

    origin = (-0.05, -0.05)
    x_end  = (0.07, -0.05)
    y_end  = (-0.05, 0.07)

    ax.annotate(
        "",
        xy=x_end, xytext=origin,
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", lw=1.8, color="black", mutation_scale=12)
    )

    ax.annotate(
        "",
        xy=y_end, xytext=origin,
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", lw=1.8, color="black", mutation_scale=12)
    )

    ax.text(
        (origin[0] + x_end[0]) / 2, origin[1] - 0.025,
        "UMAP1",
        transform=ax.transAxes,
        ha="center", va="top",
        fontsize=10
    )

    ax.text(
        origin[0] - 0.025, (origin[1] + y_end[1]) / 2,
        "UMAP2",
        transform=ax.transAxes,
        ha="right", va="center",
        rotation=90,
        fontsize=10
    )

    # category / color 가져오기
    cats = adata.obs[color].cat.categories
    colors = adata.uns[f"{color}_colors"]

    handles = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="None",
            markerfacecolor=c,
            markeredgecolor="none",
            markersize=8,
            label=lab
        )
        for lab, c in zip(cats, colors)
    ]

    leg = ax.legend(
        handles=handles,
        ncol=ncol,                     # 열 개수
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15), # 위쪽 가운데
        frameon=False,
        #fancybox=True,
        #framealpha=1,
        #facecolor="white",
        #edgecolor="white",
        fontsize=10,
        handlelength=0.8,
        handletextpad=0.4,
        columnspacing=1.4,
        borderpad=0.5
    )

    if savefig is not None:
        ax.figure.savefig(savefig, format="pdf", dpi=700, bbox_inches="tight")

    return leg


# ---------------------------------------------------------------- ORA clustermap (08_*)

def apply_style():
    """Figure defaults used by the ORA panels (Arial, editable PDF text)."""
    sns.set_style("white")
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["mathtext.fontset"] = "custom"   # 수식($...$)도 Arial로
    plt.rcParams["mathtext.rm"] = "Arial"
    plt.rcParams["mathtext.it"] = "Arial:italic"
    plt.rcParams["mathtext.bf"] = "Arial:bold"
    plt.rcParams["pdf.fonttype"] = 42             # PDF에 텍스트를 편집 가능한 형태로 임베드


def pval_to_stars(pval):
    try:
        if pval < 0.001:
            return "***"
        elif pval < 0.01:
            return "**"
        elif pval < 0.05:
            return "*"
        else:
            return ""
    except Exception:
        return ""


def row_zscore(df):
    row_mean = df.mean(axis=1)
    row_std = df.std(axis=1, ddof=0).replace(0, np.nan)
    z = df.sub(row_mean, axis=0).div(row_std, axis=0)
    return z.fillna(0)


def plot_enrichment_clustermap_combined_score_dif(
    final_results_df_dd,
    gene_sets_to_use=None,
    cluster_order=None,
    top_k=3,
    score_col="Combined Score",
    pval_col="Adjusted P-value",
    use_row_zscore=True,
    vmin=None,
    vmax=None,
    cmap="coolwarm",
    figsize=(29, 34),
    dpi=700,
    save_path=None,
    show=True,
    gene_set_palette=None,
    annot=True,
    min_score=0,
    min_specificity=None,
    tol=1e-12,
    title=None,
    title_fontsize=48,
    title_pad=30,
    cluster_color_map=None,
    col_colors_ratio=0.03,
    col_colors_gap=0.008,
    cbar_label=None,
):
    if cluster_order is None:
        cluster_order = [str(i) for i in range(17)]

    if gene_set_palette is None:
        gene_set_palette = {
            "MSigDB_Hallmark_2020": "tab:orange",
            "GO_Biological_Process_2023": "tab:blue",
            "Epithelial_state": "tab:green",
        }

    if gene_sets_to_use is None:
        gene_sets_to_use = list(gene_set_palette.keys())

    df = final_results_df_dd.copy()
    df["cluster_name"] = df["cluster_name"].astype(str)
    df = df[df["cluster_name"].isin(cluster_order)].copy()
    df = df[df["Gene_set"].isin(gene_sets_to_use)].copy()

    if df.empty:
        raise ValueError("필터링 후 남은 데이터가 없습니다. gene_sets_to_use 또는 cluster_order를 확인하세요.")

    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df[pval_col] = pd.to_numeric(df[pval_col], errors="coerce")
    df = df.dropna(subset=[score_col])

    if df.empty:
        raise ValueError(f"`{score_col}` 값이 유효한 행이 없습니다.")

    df["cluster_name"] = pd.Categorical(
        df["cluster_name"], categories=cluster_order, ordered=True
    )

    df["term_id"] = df["Gene_set"].astype(str) + " || " + df["Term"].astype(str)

    score_mat = (
        df.pivot_table(
            index="term_id",
            columns="cluster_name",
            values=score_col,
            aggfunc="max",
        )
        .reindex(columns=cluster_order)
        .fillna(0)
    )

    pval_mat = (
        df.pivot_table(
            index="term_id",
            columns="cluster_name",
            values=pval_col,
            aggfunc="min",
        )
        .reindex(columns=cluster_order)
    )

    term_meta = (
        df.sort_values(score_col, ascending=False)
        .drop_duplicates("term_id")
        .set_index("term_id")[["Term", "Gene_set"]]
    )

    selected = []

    for c in cluster_order:
        this_score = score_mat[c]

        other_cols = [x for x in cluster_order if x != c]
        if len(other_cols) > 0:
            second_score = score_mat[other_cols].max(axis=1)
        else:
            second_score = pd.Series(0, index=score_mat.index)

        specificity = this_score - second_score

        cand = pd.DataFrame({
            "owner_cluster": c,
            "term_id": score_mat.index,
            "owner_score": this_score.values,
            "second_score": second_score.values,
            "specificity": specificity.values,
        }).set_index("term_id")

        cand = cand[cand["owner_score"] > min_score].copy()
        cand = cand[cand["specificity"] >= -tol].copy()

        if min_specificity is not None:
            cand = cand[cand["specificity"] >= min_specificity].copy()

        cand = cand.sort_values(
            ["specificity", "owner_score"],
            ascending=[False, False],
        )

        cand = cand.head(top_k).reset_index()
        selected.append(cand)

    row_meta = pd.concat(selected, ignore_index=True)

    if row_meta.empty:
        raise ValueError("조건을 만족하는 top term이 없습니다.")

    row_meta["rank_in_cluster"] = row_meta.groupby("owner_cluster").cumcount() + 1
    row_meta = row_meta.join(term_meta, on="term_id")
    row_meta["row_id"] = (
        row_meta["owner_cluster"]
        + " | top"
        + row_meta["rank_in_cluster"].astype(str)
        + " | "
        + row_meta["term_id"]
    )
    row_order = row_meta["row_id"].tolist()

    raw_score_df = score_mat.reindex(
        index=row_meta["term_id"].tolist(),
        columns=cluster_order,
    ).copy()
    raw_score_df.index = row_order
    raw_score_df = raw_score_df.fillna(0)

    plot_df = raw_score_df.copy()
    if use_row_zscore:
        plot_df = row_zscore(plot_df)
        if vmin is None:
            vmin = -3
        if vmax is None:
            vmax = 3
    else:
        if vmin is None:
            vmin = 0

    annot_p = pval_mat.reindex(
        index=row_meta["term_id"].tolist(),
        columns=cluster_order,
    ).copy()
    annot_p.index = row_order

    if annot:
        annot_stars = annot_p.fillna(1.0).apply(lambda col: col.map(pval_to_stars))
    else:
        annot_stars = None

    if cbar_label is None:
        cbar_label = f" z-score ({score_col})" if use_row_zscore else score_col

    cbar_kws = {}
    if vmin is not None and vmax is not None:
        cbar_kws["ticks"] = np.linspace(vmin, vmax, 5)

    if cluster_color_map is not None:
        missing = [c for c in cluster_order if c not in cluster_color_map]
        if missing:
            raise ValueError(f"cluster_color_map에 없는 cluster가 있습니다: {missing}")

    g = sns.clustermap(
        plot_df,
        cmap=cmap,
        linecolor="white",
        linewidths=1.5,
        vmin=vmin,
        vmax=vmax,
        figsize=figsize,
        cbar_kws=cbar_kws,
        row_cluster=False,
        col_cluster=False,
        cbar_pos=(-0.03, 0.5, 0.05, 0.24),
        annot=annot_stars if annot else False,
        fmt="s",
        annot_kws={"size": 35, "color": "black"},
    )

    cbar = g.cax
    cbar.set_ylabel(cbar_label, fontsize=44)
    cbar.tick_params(labelsize=44, labelrotation=90)

    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.grid(False)

    row_to_term = row_meta.set_index("row_id")["Term"].to_dict()
    g.ax_heatmap.set_yticklabels([row_to_term[r] for r in row_order])

    # x축 label 숨김
    g.ax_heatmap.set_xticklabels([])
    g.ax_heatmap.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)

    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=37, color="black")

    if title is not None:
        g.ax_heatmap.set_title(title, fontsize=title_fontsize, pad=title_pad)

    g.ax_heatmap.margins(y=0.01)

    # heatmap 아래쪽에 annotation bar 추가
    if cluster_color_map is not None:
        heat_pos = g.ax_heatmap.get_position()
        bar_h = heat_pos.height * col_colors_ratio
        bar_y = heat_pos.y0 - col_colors_gap - bar_h

        if bar_y < 0.01:
            shift = 0.01 - bar_y
            g.ax_heatmap.set_position([
                heat_pos.x0,
                heat_pos.y0 + shift,
                heat_pos.width,
                heat_pos.height - shift,
            ])
            heat_pos = g.ax_heatmap.get_position()
            bar_h = heat_pos.height * col_colors_ratio
            bar_y = heat_pos.y0 - col_colors_gap - bar_h

        bar_colors = [cluster_color_map[c] for c in cluster_order]
        rgb = np.array([mcolors.to_rgb(c) for c in bar_colors]).reshape(1, len(cluster_order), 3)

        bar_ax = g.fig.add_axes([heat_pos.x0, bar_y, heat_pos.width, bar_h])
        bar_ax.imshow(rgb, aspect="auto", interpolation="nearest", origin="lower")
        bar_ax.set_xticks([])
        bar_ax.set_yticks([])
        for spine in bar_ax.spines.values():
            spine.set_visible(False)

    if save_path is not None:
        g.fig.savefig(save_path, format="pdf", dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(g.fig)

    return g, plot_df, row_meta, raw_score_df
