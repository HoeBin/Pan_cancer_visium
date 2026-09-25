#!/usr/bin/env python3
"""Fig 3d / Extended Data Fig 6 / Extended Data Fig 8: subtype co-enrichment networks per compartment.

Input : Source_Data edge sheets, one per network, cfg.read_source keys Fig3d (Boundary), ExtFig6_Malignant,
        ExtFig6_Normal, ExtFig8_Boundary, ExtFig8_Malignant, ExtFig8_Normal (BRCA-downsampled cohort);
        Data_Info/subtype_lineage_map.csv for node colours.
Output: Output/<Figure3|ExtendedFigure6|ExtendedFigure8>/Plots/<stem>_network_<Bdy|Mal|Normal>.pdf and
        Tables/source_data_<stem>_{edges,nodes}<suffix>.csv.
Edges are used as-is in sheet row order (the layout depends on it); nothing is recomputed. Plotting
code is that of Analysis/Run_Fig3_CoEnrichment/01_Fig3d_ExtDataFig6_network.py, so the PDFs match.
--panel fig3d|extfig6|extfig8|all.
"""
from __future__ import annotations

# "Libraries and paths" ----
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # Code/ for `common`
from common import config as cfg

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.colors import Normalize

plt.rcParams.update({
    "font.family": "Arial", "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})

# "Inputs" ----
LINEAGE_COLORS = cfg.LINEAGE_COLORS

# Manuscript edge rule; used only for the title and the colour-scale floor (the edges come from the sheet).
ALPHA, VOTE, K_MIN, REP_Q, J_MIN = 0.01, "p", 5, 0.05, 0.5001

# --panel -> [(sheet key, compartment code, output figure, file stem, table suffix)]
# files: <stem>_network_<code>.pdf, source_data_<stem>_{edges,nodes}<suffix>.csv
PANELS = {
    "fig3d": [("Fig3d", "Bdy", "Fig3", "fig3d", "")],
    "extfig6": [("ExtFig6_Malignant", "Mal", "ExtFig6", "extfig6", "_Malignant"),
                ("ExtFig6_Normal", "Normal", "ExtFig6", "extfig6", "_Normal")],
    "extfig8": [("ExtFig8_Boundary", "Bdy", "ExtFig8", "extfig8", "_Boundary"),
                ("ExtFig8_Malignant", "Mal", "ExtFig8", "extfig8", "_Malignant"),
                ("ExtFig8_Normal", "Normal", "ExtFig8", "extfig8", "_Normal")],
}


def load_lineage_map() -> dict[str, str]:
    """subtype -> major cell-type group (one of the 7 manuscript lineages)."""
    return cfg.load_subtype_lineage()                       # Data_Info/subtype_lineage_map.csv


# "Layout and drawing (same code as Analysis/Run_Fig3_CoEnrichment/01_Fig3d_ExtDataFig6_network.py)" ----
def ordered_subgraph(G: nx.Graph, nodes: list) -> nx.Graph:
    """Induced subgraph with node order exactly `nodes`; nx.subgraph() may iterate its node set
    (hash order), which would make the layout process-dependent."""
    H = nx.Graph()
    H.add_nodes_from(nodes)
    keep = set(nodes)
    H.add_edges_from((u, v, d) for u, v, d in G.edges(data=True) if u in keep and v in keep)
    return H


def fig3d_layout(G: nx.Graph, node_lineage: dict[str, str], seed: int) -> dict:
    """Two-block layout: TME band on top, epithelial block bottom-right, gateway subtypes
    (TME nodes with an epithelial edge) on the corridor between the two blocks."""
    epi = [n for n in G if node_lineage[n] == "Epithelial"]
    tme = [n for n in G if node_lineage[n] != "Epithelial"]
    if not epi or not tme:
        return nx.spring_layout(G, weight="weight", seed=seed,
                                k=1.6 / max(1, np.sqrt(G.number_of_nodes())))

    def clamp(pos: dict, q: float = 0.85) -> dict:
        """Pull outliers beyond the q-quantile radius back onto that radius."""
        P = np.array(list(pos.values()), float)
        c = P.mean(axis=0)
        d = np.linalg.norm(P - c, axis=1)
        r = float(np.quantile(d, q)) or 1.0
        out = {}
        for n, p in pos.items():
            v = np.asarray(p, float) - c
            dist = np.linalg.norm(v) or 1e-9
            out[n] = tuple(c + v * min(1.0, r / dist))
        return out

    def rescale(pos: dict, x0: float, x1: float, y0: float, y1: float) -> dict:
        xs = np.array([p[0] for p in pos.values()])
        ys = np.array([p[1] for p in pos.values()])
        xr = (xs.max() - xs.min()) or 1.0
        yr = (ys.max() - ys.min()) or 1.0
        return {n: ((p[0] - xs.min()) / xr * (x1 - x0) + x0,
                    (p[1] - ys.min()) / yr * (y1 - y0) + y0) for n, p in pos.items()}

    def pca_align(pos: dict) -> dict:
        """Rotate so the block's principal axis becomes horizontal (SVD sign
        fixed deterministically so reruns cannot mirror the block)."""
        P = np.array(list(pos.values()), float)
        X = P - P.mean(axis=0)
        _, _, vt = np.linalg.svd(X, full_matrices=False)
        Y = X @ vt.T
        for k in range(Y.shape[1]):
            j = int(np.argmax(np.abs(Y[:, k])))
            if Y[j, k] < 0:
                Y[:, k] *= -1
        return {n: tuple(y) for n, y in zip(pos, Y)}

    bridges = sorted(n for n in tme if any(node_lineage[m] == "Epithelial" for m in G[n]))
    pos_t = nx.spring_layout(ordered_subgraph(G, tme), weight="weight", seed=seed,
                             k=1.5 / np.sqrt(len(tme)))
    pos_e = nx.spring_layout(ordered_subgraph(G, epi), weight="weight", seed=seed,
                             k=1.7 / np.sqrt(max(2, len(epi))))

    # Γ-shaped composition: TME band on top, epithelial block below its right arm
    pos_t = pca_align(clamp(pos_t))
    imm = [n for n in tme if node_lineage[n] in ("T cell", "B cell", "Myeloid")]
    str_ = [n for n in tme if node_lineage[n] in ("Fibroblast", "Endothelial", "Mural")]
    if imm and str_:
        if (np.mean([pos_t[n][0] for n in imm]) < np.mean([pos_t[n][0] for n in str_])):
            pos_t = {n: (-p[0], p[1]) for n, p in pos_t.items()}
    pos = rescale(pos_t, -1.05, 1.05, 0.10, 1.05)
    pos.update(rescale(pca_align(clamp(pos_e)), -0.10, 1.10, -1.75, -0.95))

    # gateways on the corridor between the blocks, spread sideways
    c_src = np.mean([pos[n] for n in (imm or tme) if n not in bridges], axis=0)
    c_e = np.mean([pos[n] for n in epi], axis=0)
    d = c_e - c_src
    perp = np.array([-d[1], d[0]]) / (np.linalg.norm(d) or 1.0)
    for i, n in enumerate(bridges):
        t = 0.45 + 0.28 * (i + 1) / (len(bridges) + 1)
        pos[n] = tuple(c_src + t * d + perp * 0.22 * (1 if i % 2 == 0 else -1))
    # cluster each gateway's strongest partners around it so bridges stay short
    placed = set(bridges)
    for n in bridges:
        strong = sorted((m for m in G[n] if node_lineage[m] != "Epithelial"
                         and m not in bridges and G[n][m]["weight"] >= 0.535
                         and G.degree(m) <= 20),
                        key=lambda m: -G[n][m]["weight"])
        for j, m in enumerate(strong):
            side = 0 if j == 0 else (j + 1) // 2 * (1 if j % 2 == 0 else -1)
            off = np.array([0.36 * side, 0.30 + 0.12 * (j % 2)])
            pos[m] = tuple(np.asarray(pos[n]) + off)
            placed.add(m)
        eps = sorted(m for m in G[n] if node_lineage[m] == "Epithelial")
        for j, m in enumerate(eps):
            off = np.array([0.20 * (j - (len(eps) - 1) / 2), -0.24])
            pos[m] = tuple(np.asarray(pos[n]) + off)
            placed.add(m)

    # spread the dense immune cluster vertically so labels clear
    core_imm = [n for n in imm if pos[n][1] > 0.05]
    if core_imm:
        cy = float(np.mean([pos[n][1] for n in core_imm]))
        for n in core_imm:
            x, y = pos[n]
            pos[n] = (x, cy + (y - cy) * 1.45)

    # pull sparse peripheral nodes toward their neighbors' centroid
    for n in G.nodes:
        if n in placed or node_lineage[n] == "Epithelial" or G.degree(n) > 5:
            continue
        nb = [pos[m] for m in G[n]]
        if nb:
            c = np.mean(nb, axis=0)
            pos[n] = tuple(0.3 * np.asarray(pos[n]) + 0.7 * c)
    return pos


def draw_spread_labels(ax, pos: dict, fontsize: float, fig_w_in: float, iters: int = 300) -> None:
    """Node labels anchored at node centers; overlapping label boxes are pushed
    apart pairwise (with a weak spring back to the node) until they clear."""
    names = list(pos)
    P = np.array([pos[n] for n in names], float)
    span = max(np.ptp(P[:, 0]), np.ptp(P[:, 1])) or 1.0
    pt = span / (fig_w_in * 72 * 0.75)  # approx. data units per point
    hw = np.array([0.5 * (len(n) + 1) * fontsize * 0.56 * pt for n in names])
    hh = np.full(len(names), 0.62 * fontsize * pt)
    L = P.copy()
    for _ in range(iters):
        moved = False
        for i in range(len(names) - 1):
            for j in range(i + 1, len(names)):
                dx, dy = L[j] - L[i]
                ox = hw[i] + hw[j] - abs(dx)
                oy = hh[i] + hh[j] - abs(dy)
                if ox > 0 and oy > 0:
                    moved = True
                    if ox * 0.3 < oy:
                        s = 0.5 * ox * (1 if dx >= 0 else -1)
                        L[j, 0] += s
                        L[i, 0] -= s
                    else:
                        s = 0.5 * oy * (1 if dy >= 0 else -1)
                        L[j, 1] += s
                        L[i, 1] -= s
        # the node center may not leave its label box
        D = L - P
        D[:, 0] = np.clip(D[:, 0], -hw * 0.9, hw * 0.9)
        D[:, 1] = np.clip(D[:, 1], -hh * 2.2, hh * 2.2)
        L = P + D
        L += (P - L) * 0.02
        if not moved:
            break
    for k, n in enumerate(names):
        ax.text(L[k, 0], L[k, 1], n, fontsize=fontsize, ha="center", va="center",
                color="#111111", zorder=5)


def draw_network(sheet_key: str, location: str, figure: str, stem: str, sfx: str, seed: int) -> None:
    """One network: read the edge sheet, build the graph in sheet row order, write the edge/node CSVs, draw the PDF."""
    figdir, tabdir = cfg.output_dirs(figure)
    edges = cfg.read_source(sheet_key)                  # sheet rows as-is
    print(f"[{sheet_key}] edges={len(edges)}  (k>={K_MIN}, J>={J_MIN}, rep_q<{REP_Q}; taken from the sheet)")

    lineage_map = load_lineage_map()

    # edge (and node) insertion order = sheet row order; the layout depends on it
    G = nx.Graph()
    for _, e in edges.iterrows():
        # 15 significant digits (xlsx precision) so both layers get identical layout input
        G.add_edge(e.node_a, e.node_b, weight=cfg.excel15(e.pan_cancer_median_jaccard),
                   n_sig=int(e.n_cancers_significant))
    print(f"nodes={G.number_of_nodes()}  isolated subtypes excluded={93 - G.number_of_nodes()}")

    node_lineage = {n: lineage_map.get(n, "Myeloid") for n in G.nodes}
    nodes_df = pd.DataFrame({
        "node": list(G.nodes),
        "lineage": [node_lineage[n] for n in G.nodes],
        "degree": [G.degree(n) for n in G.nodes],
    })
    edges.to_csv(tabdir / f"source_data_{stem}_edges{sfx}.csv", index=False)
    nodes_df.to_csv(tabdir / f"source_data_{stem}_nodes{sfx}.csv", index=False)

    pos = fig3d_layout(G, node_lineage, seed)

    fig, ax = plt.subplots(figsize=(12.5, 13.5))
    wts = np.array([G[u][v]["weight"] for u, v in G.edges])
    norm = Normalize(vmin=max(J_MIN, wts.min()), vmax=np.quantile(wts, 0.98))
    cmap = (matplotlib.colormaps["YlOrRd"] if hasattr(matplotlib, "colormaps") else cm.get_cmap("YlOrRd"))
    ecolors = [cmap(0.35 + 0.65 * norm(w)) for w in wts]
    nx.draw_networkx_edges(G, pos, ax=ax, width=0.6 + 2.4 * norm(wts), edge_color=ecolors, alpha=0.75)
    deg = np.array([G.degree(n) for n in G.nodes], dtype=float)
    sizes = 60 + 240 * (deg / max(1.0, deg.max()))
    ncolors = [LINEAGE_COLORS[node_lineage[n]] for n in G.nodes]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=sizes, node_color=ncolors,
                           edgecolors="#333333", linewidths=0.5)
    draw_spread_labels(ax, pos, fontsize=13.0, fig_w_in=12.5)

    vote_txt = f"sign-test p<{ALPHA}" if VOTE == "p" else f"BH q<{ALPHA}"
    ax.set_title(
        f"{location} co-enrichment network — spatial-null, replicated "
        f"(vote: {vote_txt}; rep-BH q<{REP_Q}; ≥{K_MIN} cancers; pan-median J≥{J_MIN})",
        fontsize=13,
    )
    ax.axis("off")
    fig.tight_layout()
    final_pdf = figdir / f"{stem}_network_{location}.pdf"
    fig.savefig(final_pdf, bbox_inches="tight")
    plt.close(fig)
    print("saved:", final_pdf)


def main() -> None:
    ap = argparse.ArgumentParser(description="Fig 3d / Extended Data Fig 6 / Extended Data Fig 8 co-enrichment "
                                             "networks from the Source Data edge sheets.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--panel", choices=("fig3d", "extfig6", "extfig8", "all"), default="all",
                    help="which networks to draw (fig3d: Boundary; extfig6: Malignant/Normal; "
                         "extfig8: BRCA-downsampled Boundary/Malignant/Normal)")
    ap.add_argument("--seed", type=int, default=42, help="spring_layout seed")
    args = ap.parse_args()

    panels = list(PANELS) if args.panel == "all" else [args.panel]
    for p in panels:
        for sheet_key, location, figure, stem, sfx in PANELS[p]:
            draw_network(sheet_key, location, figure, stem, sfx, args.seed)


if __name__ == "__main__":
    main()
