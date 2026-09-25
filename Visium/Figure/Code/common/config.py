# Shared settings for the Python Figure scripts (Fig.3/4/5, ExtFig.3): paths, palettes, orders, Source Data
# readers. Counterpart of common/config.R for the R scripts. Scripts import it with
#   sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
#   from common import config as cfg
# Figure-layer only: this repository does not include the Analysis layer (raw h5ad / TCGA counts / bundled
# statistics Input) of the source project. The only bundled non-Source-Data inputs kept here are the two small
# files the h5ad-dependent image panels need (Figure4/Input/pemt_scores.csv.gz, Figure5/Input/selected_slides.csv).

from __future__ import annotations

import os
from pathlib import Path

# "Paths" ---------------------------------------------------------------------
CODE_DIR = Path(__file__).resolve().parents[1]        # .../Github_code_v2/Code
ROOT = CODE_DIR.parent                                  # .../Github_code_v2
DATA_INFO_DIR = CODE_DIR / "Data_Info"
SOURCE_DATA_DIR = ROOT / "Paper_info" / "Paper_info" / "Source_Data"   # same files common/config.R reads
OUTPUT_DIR = ROOT / "Output"                                            # same root the R scripts write to

SUBTYPE_LINEAGE_CSV = DATA_INFO_DIR / "subtype_lineage_map.csv"  # 93 subtypes -> 7 lineages (Fig3d/ExtFig6/8 node colours)

# Small bundled files the h5ad-dependent panels need (not Source Data, kept alongside their scripts).
FIG4E_INPUT = CODE_DIR / "Figure4" / "Input"   # pemt_scores.csv.gz
FIG5C_INPUT = CODE_DIR / "Figure5" / "Input"   # selected_slides.csv


def _env_path(var: str, default: str) -> Path:
    """Path from the environment variable var, or the placeholder default."""
    return Path(os.environ.get(var, default))


# Raw data (not bundled): the deposited Visium AnnData, needed only by the three h5ad image panels
# (Figure4/06, Figure4/07, Figure5/03). Set PCASSO_VISIUM_H5AD or pass --h5ad / --h5ad-path.
VISIUM_H5AD = _env_path("PCASSO_VISIUM_H5AD", "/path/to/_visium_ALL_260730.h5ad")


def require_external(path: Path, what: str) -> Path:
    """Check that an unbundled raw file exists; exit with instructions otherwise."""
    if str(path).startswith("/path/to/") or not Path(path).exists():
        raise SystemExit(
            f"[config] {what} not found: {path}\n"
            f"         This file is not bundled. Set PCASSO_VISIUM_H5AD or pass the path on the command line."
        )
    return Path(path)


# "Output directories" --------------------------------------------------------
FIGURE_DIRS = {
    "Fig3": "Figure3", "Fig4": "Figure4", "Fig5": "Figure5",
    "ExtFig3": "ExtendedFigure3", "ExtFig5": "ExtendedFigure5", "ExtFig6": "ExtendedFigure6", "ExtFig7": "ExtendedFigure7",
    "ExtFig8": "ExtendedFigure8", "ExtFig9": "ExtendedFigure9", "ExtFig10": "ExtendedFigure10",
    "SuppFig11": "SupplementaryFigure11", "SuppFig12_13": "SupplementaryFigure12_13",
}


def output_dirs(figure: str) -> tuple[Path, Path]:
    """Create Output/<figure>/Plots and Output/<figure>/Tables; return (plots, tables)."""
    base = OUTPUT_DIR / FIGURE_DIRS.get(figure, figure)
    plots, tables = base / "Plots", base / "Tables"
    plots.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    return plots, tables


def intermediate_dir(figure: str, *sub: str) -> Path:
    """Create and return Output/<figure>/Intermediate[/sub...] (files re-read by other scripts)."""
    d = OUTPUT_DIR / FIGURE_DIRS.get(figure, figure) / "Intermediate"
    for s in sub:
        d = d / s
    d.mkdir(parents=True, exist_ok=True)
    return d


# "Source Data (Figure-layer input)" -------------------------------------------
# Counterpart of read_source(fig_file, sheet) in common/config.R.
SOURCE_FILES = {
    "Fig1": "Source_Data_Fig1.xlsx", "Fig2": "Source_Data_Fig2.xlsx", "Fig3": "Source_Data_Fig3.xlsx",
    "Fig4": "Source_Data_Fig4.xlsx", "Fig5": "Source_Data_Fig5.xlsx",
    "ExtFig1": "Source_Data_Extended_Fig1.xlsx", "ExtFig2": "Source_Data_Extended_Fig2.xlsx",
    "ExtFig3": "Source_Data_Extended_Fig3.xlsx", "ExtFig4": "Source_Data_Extended_Fig4.xlsx",
    "ExtFig5": "Source_Data_Extended_Fig5.xlsx", "ExtFig6": "Source_Data_Extended_Fig6.xlsx",
    "ExtFig7": "Source_Data_Extended_Fig7.xlsx", "ExtFig8": "Source_Data_Extended_Fig8.xlsx",
    "ExtFig10": "Source_Data_Extended_Fig10.xlsx", "Supp": "Source_Data_Supp_Fig.xlsx",
}
# Panel key -> (file key, sheet name, header row). header=1: the first row is a title (e.g. "Malignant");
# header=None: several titled sections in one sheet, read with read_source_sections().
SHEETS = {
    "Fig3b": ("Fig3", "Fig.3b", 0), "Fig3c": ("Fig3", "Fig.3c", 0), "Fig3d": ("Fig3", "Fig.3d", 0),
    "Fig3e": ("Fig3", "Fig.3e", 0), "Fig3f": ("Fig3", "Fig.3f", 0), "Fig3g": ("Fig3", "Fig.3g", 0),
    "Fig3h": ("Fig3", "Fig.3h", 0),
    "Fig4a": ("Fig4", "Fig.4a", 0), "Fig4b": ("Fig4", "Fig.4b", 0), "Fig4d": ("Fig4", "Fig.4d", 0),
    "Fig4f": ("Fig4", "Fig.4f", 0), "Fig4g": ("Fig4", "Fig.4g", 0),
    "Fig5d_left": ("Fig5", "Fig.5d left", 0), "Fig5d_right": ("Fig5", "Fig.5d right", 0),
    "ExtFig5b": ("ExtFig5", "ED Fig.5b", 0),
    "ExtFig6_Malignant": ("ExtFig6", "ED Fig.6 top", 1), "ExtFig6_Normal": ("ExtFig6", "ED Fig.6 bottom", 1),
    "ExtFig7b": ("ExtFig7", "ED Fig.7b", 0), "ExtFig7c": ("ExtFig7", "ED Fig.7c", 0),
    "ExtFig7d": ("ExtFig7", "ED Fig.7d", 0), "ExtFig7e": ("ExtFig7", "ED Fig.7e", 0),
    "ExtFig7f": ("ExtFig7", "ED Fig.7f", 0), "ExtFig7g": ("ExtFig7", "ED Fig.7g", 0),
    "ExtFig7h": ("ExtFig7", "ED Fig.7h", 0),
    "ExtFig8_Malignant": ("ExtFig8", "ED Fig.8 top", 1), "ExtFig8_Boundary": ("ExtFig8", "ED Fig.8 middle", 1),
    "ExtFig8_Normal": ("ExtFig8", "ED Fig.8 bottom", 1),
    "ExtFig3": ("ExtFig3", "ED Fig.3", 0),
    "ExtFig10": ("ExtFig10", "ED Fig.10", 0),
    "SuppFig11": ("Supp", "Sup fig.11", 0),
    "SuppFig12_top": ("Supp", "Sup fig.12 top", None), "SuppFig12_bottom": ("Supp", "Sup fig. 12 bottom", None),
    "SuppFig13": ("Supp", "Sup fig.13", None),
}

# Aliases to the Source Data spelling (subtype `Fibroblast`, γδ T as `γδ_T`).
SUBTYPE_ALIASES = {"gd_T": "γδ_T", "Normal_Fibroblast": "Fibroblast", "Normal Fibroblast": "Fibroblast",
                   "款灌 T": "γδ T"}   # ED Fig.3 sheet spelling of "γδ T" (UTF-8 bytes read as CP949); kept as a safeguard


def excel15(x):
    """Truncate (not round) a float to 15 significant digits, the precision Excel writes to xlsx.
    Source Data value == excel15(pipeline value), so value-sensitive steps (e.g. spring_layout) use it in both layers."""
    import math
    from decimal import Decimal, ROUND_DOWN
    x = float(x)
    if x == 0.0 or not math.isfinite(x):
        return x
    d = Decimal(repr(x))
    return float(d.quantize(Decimal(1).scaleb(d.adjusted() - 14), rounding=ROUND_DOWN))


def canonical_subtype(name):
    """Map a subtype name to its Source Data spelling (unchanged if not an alias)."""
    return SUBTYPE_ALIASES.get(name, name) if isinstance(name, str) else name


def canonical_pair(pair: str, sep: str = "@") -> str:
    """Apply canonical_subtype to both parts of an 'a@b' subtype-pair string."""
    return sep.join(canonical_subtype(x) for x in pair.split(sep))


def source_path(fig_file: str) -> Path:
    """File key ('Fig3') or file name ('Source_Data_Fig3.xlsx') -> path under Paper_info/Paper_info/Source_Data/."""
    name = SOURCE_FILES.get(fig_file, fig_file)
    path = SOURCE_DATA_DIR / name
    if not path.exists():
        raise SystemExit(f"[config] Source Data file not found: {path}")
    return path


def _tidy(df):
    df = df.loc[:, [c for c in df.columns if not str(c).startswith("Unnamed")]]
    df = df.dropna(how="all").reset_index(drop=True)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def read_source(key_or_file: str, sheet: str | None = None, header: int = 0):
    """Read one Source Data sheet: read_source("Fig3d") (panel key), read_source("Fig3", "Fig.3d") or
    read_source("Source_Data_Fig3.xlsx", "Fig.3d"). Drops Unnamed columns and empty rows; strips column names."""
    import pandas as pd
    if sheet is None:
        if key_or_file not in SHEETS:
            raise KeyError(f"unknown Source Data panel key: {key_or_file!r}")
        fig_file, sheet, header = SHEETS[key_or_file]
        if header is None:
            raise ValueError(f"{key_or_file!r} has several sections; use read_source_sections()")
    else:
        fig_file = key_or_file
    df = pd.read_excel(source_path(fig_file), sheet_name=sheet, header=header, engine="openpyxl")
    return _tidy(df)


def read_source_sections(key_or_file: str, sheet: str | None = None) -> dict:
    """Read a sheet split into titled sections (e.g. Sup fig.12: 'Cox proportional hazard model', 'Kaplan-meier')
    as {title: DataFrame}. A title row has only its first cell filled; the next row is the table header."""
    import openpyxl
    import pandas as pd
    if sheet is None:
        fig_file, sheet, _ = SHEETS[key_or_file]
    else:
        fig_file = key_or_file
    ws = openpyxl.load_workbook(source_path(fig_file), read_only=True)[sheet]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    out, i = {}, 0
    while i < len(rows):
        r = rows[i]
        nonempty = [v for v in r if v is not None]
        if len(nonempty) == 1 and r[0] is not None and i + 1 < len(rows):
            title = str(r[0]).strip()
            hdr = [str(v).strip() for v in rows[i + 1] if v is not None]
            data = []
            j = i + 2
            while j < len(rows) and any(v is not None for v in rows[j]):
                data.append(rows[j][:len(hdr)])
                j += 1
            out[title] = _tidy(pd.DataFrame(data, columns=hdr))
            i = j
        else:
            i += 1
    return out


# "Loaders" -------------------------------------------------------------------
def load_subtype_lineage() -> dict[str, str]:
    """subtype_lineage_map.csv -> {subtype: lineage} (93 subtypes -> 7 lineages)."""
    import pandas as pd
    mp = pd.read_csv(SUBTYPE_LINEAGE_CSV, dtype=str)
    return dict(zip(mp["subtype"], mp["lineage"]))


# "Palettes and orders" -------------------------------------------------------
# Myeloid/Fibroblast/Mural colors differ from the reference repository's config.R.
LINEAGE_ORDER = ["Epithelial", "T cell", "B cell", "Myeloid", "Endothelial", "Fibroblast", "Mural"]
LINEAGE_COLORS = {
    "Epithelial": "#FF0000", "T cell": "#0000FF", "B cell": "#00FF00", "Myeloid": "#9D00FF",
    "Endothelial": "#FF00FF", "Fibroblast": "#FF9500", "Mural": "#00FFFF",
}

COMP_ORDER = ["Mal", "Bdy", "Normal"]                                   # compartment codes (obs LocationNew)
COMP_LABELS = {"Mal": "Malignant", "Bdy": "Boundary", "Normal": "Normal"}
COMP_COLORS = {"Mal": "#D62728", "Bdy": "#2CA02C", "Normal": "#4477CC"}
PAIR_CATEGORIES = ["Epi-Epi", "Epi-TME", "TME-TME"]                     # subtype-pair category (ExtFig5b, ExtFig7d)

CANCER_ORDER = ["BRCA", "COCA", "HNCA", "KICA", "LICA", "LUCA", "OVCA", "PACA", "PRCA", "SKCA", "STCA", "THCA", "UECA"]
# Fig3b per-cancer (marker, color, open); open=True draws a hollow marker.
CANCER_STYLE = {
    "BRCA": ("o", "#7570B3", True), "COCA": ("^", "#E7298A", True), "HNCA": ("+", "#A6761D", False),
    "KICA": ("x", "#666666", False), "LICA": ("D", "#5C4D8F", True), "LUCA": ("v", "#FF57C0", True),
    "OVCA": ("s", "#8DC244", True), "PACA": ("*", "#F5C832", False), "PRCA": ("d", "#8A5510", True),
    "SKCA": ("P", "#8A8A8A", False), "STCA": ("h", "#A999D7", True), "THCA": ("p", "#CDE67F", True),
    "UECA": ("X", "#66A61E", False),
}


# "Extended Data Fig. 3 (marker gene x deconvolution correlation heatmaps)" ---------
# Per panel: (display name, obs deconvolution-abundance column, marker gene). List order = heatmap column order;
# unique genes in order = row order (a gene shared by several subtypes gets one row, e.g. SPP1).
EXTFIG3_PANEL_LETTER = {"Major celltypes": "a", "B cell": "b", "Myeloid": "c", "T_NK cell": "d",
                        "Endothelial": "e", "Fibroblast": "f", "Epithelial": "g", "Mural": "h"}
EXTFIG3_MAJORS = [("Epithelial", "CDH1"), ("T cell", "CD3E"), ("B cell", "CD79A"), ("Myeloid", "CD68"),
                  ("Endothelial", "RAMP2"), ("Fibroblast", "COL1A1"), ("Mural", "ACTA2")]   # (major, marker); obs column = "<major>_enriched"
EXTFIG3_PANELS = [
    ("T_NK cell", [
        ("CD8 Tcm/Tn", "CD8 Tcm/Tn", "TCF7"), ("CD8 Tem/Trm", "CD8 Tem/Trm", "SH2D1A"),
        ("CD8 Tem/Temra", "CD8 Tem/Temra", "GZMH"), ("CD8 Trm", "CD8 Trm", "ITGA1"),
        ("CD8 Tex", "CD8 Tex", "LAG3"), ("CD8 Cycling Tex", "CD8 Cycling Tex", "HAVCR2"),
        ("ISG15+ CD8 T", "CD8 ISG15+ T", "IFI44"), ("CD4 Tcm/Tn", "CD4 Tcm/Tn", "LEF1"),
        ("CD4 Tem/Effector", "CD4 Tem/Effector", "AQP3"), ("CD4 Th1", "CD4 Th1", "ANXA1"),
        ("CD4 Th17", "CD4 Th17", "IL7R"), ("CD4 Tfh", "CD4 Tfh", "AREG"), ("CD4 Treg", "CD4 Treg", "CD4"),
        ("CD4 Cycling Treg", "CD4 Cycling Treg", "FOXP3"), ("CD4 CXCL13+ T", "CD4 CXCL13+ T", "CXCL13"),
        ("NK", "NK", "NKG7"), ("CD16+ NK", "CD16+ NK", "IGFBP7"), ("CD16- NK", "CD16- NK", "IRF8"),
        ("γδ T", "γδ T", "KIR2DL4"), ("MAIT", "MAIT", "SLC4A10"), ("ILC3", "ILC3", "KIT"),
    ]),
    ("B cell", [
        ("Bn", "Bn", "CD72"), ("NR4A2+ Bn", "NR4A2+ Bn", "CREM"), ("Bgc", "Bgc", "CD40"),
        ("Bmem", "Bmem", "GPR183"), ("ABC", "ABC", "HCK"), ("Plasma", "Plasma", "MZB1"),
        ("Cycling Bgc", "Cycling Bgc", "TOP2A"), ("Plasmablast", "Plasmablast", "SDC1"),
    ]),
    ("Myeloid", [
        ("CD14+Mono", "CD14+_Mono", "S100A12"), ("CD16+Mono", "CD16+_Mono", "LILRB2"), ("Macro", "Macro", "C1QC"),
        ("CXCL3+ Macro", "CXCL3+_Macro", "CXCL3"), ("FOLR2+ Macro", "FOLR2+_Macro", "FOLR2"),
        ("HSPA6+ Macro", "HSPA6+_Macro", "HSPA6"), ("SPP1+ Macro", "SPP1+_Macro", "SPP1"),
        ("SPP1+ CXCL3+ Macro", "SPP1+_CXCL3+_Macro", "SPP1"), ("SPP1+ MT1+ Macro", "SPP1+_MT1+_Macro", "MT1F"),
        ("SPP1+ cycMacro", "SPP1+_cycMacro", "MKI67"), ("TNFSF10+ Macro", "TNFSF10+_Macro", "TNFSF10"),
        ("LAMP3+ cDC", "LAMP3+_cDC", "LAMP3"), ("Mast", "Mast", "KIT"),
    ]),
    ("Fibroblast", [
        ("apCAF", "apCAF", "CD74"), ("iCAF", "iCAF", "C7"), ("IL6+ iCAF", "IL6+_iCAF", "IL6"),
        ("Pl16+ iCAF", "PI16+_iCAF", "CFD"), ("myoCAF", "myoCAF", "PLN"), ("pnCAF", "pnCAF", "GPM6B"),
        ("vCAF", "vCAF", "COL18A1"), ("mCAF", "mCAF", "POSTN"), ("cyc mCAF", "cyc_mCAF", "TOP2A"),
        ("CXCL14+ mCAF", "CXCL14+_mCAF", "CXCL14"), ("ISG15 + mCAF", "ISG15+_mCAF", "ISG15"),
        ("tCAF", "tCAF", "NDRG1"), ("HSP+ tCAF", "HSP+_tCAF", "PDPN"),
    ]),
    ("Endothelial", [
        ("Arterial EC", "Arterial_EC", "GJA5"), ("Capillary EC", "Capillary_EC", "RGCC"),
        ("Venous EC", "Venous_EC", "CLU"), ("Venous iEC", "Venous_iEC", "SELE"),
        ("Lymphatic EC", "Lymphatic_EC", "PDPN"), ("Tip EC", "Tip_EC", "KDR"),
        ("PGF+ Tip EC", "PGF+_Tip_EC", "COL4A1"), ("Fibrosis PGF+ TipEC", "Fibrosis_PGF+_Tip_EC", "NOTCH3"),
        ("TMEM100+ EC", "TMEM100+_EC", "TMEM100"), ("SPRY1+ EC", "SPRY1+_EC", "SPRY1"),
        ("HMOX1+ iEC", "HMOX1+_iEC", "HMOX1"), ("CCL2+ iEC", "CCL2+_iEC", "CCL2"),
        ("CD14+ EC", "CD14+ EC", "FCN3"), ("HMGB2+ cycEC", "HMGB2+_cycEC", "TOP2A"),
    ]),
    ("Epithelial", [
        ("Epi_Cycle", "Epi_Cycle", "TOP2A"), ("Epi_Luminal", "Epi_Luminal", "ESR1"),
        ("Epi_Hormonal", "Epi_Hormonal", "MLPH"), ("Epi_Xenobiotic", "Epi_Xenobiotic", "AKR1C2"),
        ("Epi_TNFinfla", "Epi_TNF", "TNFRSF12A"), ("Epi_Interferon", "Epi_Interferon", "IFIT3"),
        ("Epi_Glandular", "Epi_Glandular", "ELF3"), ("Epi_Squamous", "Epi_Squamous", "CSTA"),
        ("Epi_Hypox-Stress", "Epi_Hypox-Stress", "NUPR1"), ("Epi_cEMT", "Epi_cEMT", "MARCKSL1"),
        ("Epi_ER-Stress", "Epi_ER-Stress", "CLU"), ("Epi_Basal", "Epi_Basal", "CALD1"),
        ("Epi_Hypox-Adapt", "Epi_Hypox-Adapt", "NDUFA4L2"), ("Epi_Ciliated", "Epi_Ciliated", "TPPP3"),
        ("Epi_Oxphos-Ion", "Epi_Oxphos-Ion", "SCIN"), ("Epi_Oxphos-Metal", "Epi_Oxphos-Metal", "GPX3"),
    ]),
    ("Mural", [("SMC", "SMC", "MYL9"), ("Pericyte", "Pericyte", "RGS5")]),
]

# "matplotlib defaults" -------------------------------------------------------
def setup_matplotlib(font: str = "Arial"):
    """Agg backend, Arial (or the default sans-serif) font, editable text (pdf/ps fonttype 42)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    try:
        font_manager.findfont(font, fallback_to_default=False)
        family = font
    except ValueError:
        family = plt.rcParams["font.sans-serif"][0]
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [family, "Liberation Sans", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    return family
