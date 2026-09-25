#!/usr/bin/env python3
"""Ext. Data Fig. 9: the nine per-slide PDFs of 05_Fig4e_ExtDataFig9_slide_panels.py
tiled on one page, 3 x 3, unscaled.

Input : Output/Figure4/Plots/Fig4e_ExtDataFig9_slide_panels/<cancer_type>/<sample_id>__pemt_overlap_panels.pdf
        for DEFAULT_SAMPLES (row by row), or rows x cols single-page PDFs given as
        positional arguments. Run 05 first.
Output: Output/ExtendedFigure9/Plots/ExtDataFig9_grid_3x3_9cancers.pdf (--out-pdf).
Pages are copied with pypdf (vector spots, H&E images and fonts unchanged); every cell takes the largest page
width and height so pages never overlap, each page top-aligned (--align) and centred in its cell.
Options: --rows/--cols (3 x 3), --gap and --margin in inches (0.5 each). Needs pypdf.
"""
from __future__ import annotations

# "Libraries and paths" ----
import argparse
import sys
from pathlib import Path
from typing import List

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # Code/ for `common`
from common import config as cfg

from pypdf import PageObject, PdfReader, PdfWriter, Transformation

PT_PER_INCH = 72.0

# "Inputs" ----
# Per-slide PDFs written by 05_Fig4e_ExtDataFig9_slide_panels.py (its default --outdir).
PANEL_DIR = cfg.output_dirs("Fig4")[0] / "Fig4e_ExtDataFig9_slide_panels"

# "Outputs" ----
DEFAULT_OUT_PDF = cfg.output_dirs("ExtFig9")[0] / "ExtDataFig9_grid_3x3_9cancers.pdf"

# The nine slides of Extended Data Fig. 9, in grid order (row by row).
DEFAULT_SAMPLES = [
    "scCRLM_Atlas_ST-colon1",
    "GSM5924042_frozen_a_1",
    "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_220714",
    "GSM6506110_SP1",
    "GSM7498812_SS1923404",
    "GSM7211257_EHU-W3",
    "GSE251950_20_00331_LI_SING",
    "S20_63981",
    "01_034_C3d1",
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Tile single-page PDFs on one page in a rows x cols grid, unscaled.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("pdfs", nargs="*", help="Input single-page PDFs, in grid order (row by row). "
                        "Default: the nine Extended Data Fig. 9 panels under "
                        "Output/Figure4/Plots/Fig4e_ExtDataFig9_slide_panels/.")
    parser.add_argument("--out-pdf", default=str(DEFAULT_OUT_PDF),
                        help="Output PDF path.")
    parser.add_argument("--rows", type=int, default=3)
    parser.add_argument("--cols", type=int, default=3)
    parser.add_argument("--gap", type=float, default=0.5, help="Gap between cells (inches).")
    parser.add_argument("--margin", type=float, default=0.5, help="Page margin (inches).")
    parser.add_argument("--align", choices=["top", "center"], default="top",
                        help="Vertical alignment of each page inside its cell.")
    args = parser.parse_args()

    if not args.pdfs:
        args.pdfs = []
        for sample in DEFAULT_SAMPLES:
            hits = sorted(PANEL_DIR.glob(f"*/{sample}__pemt_overlap_panels.pdf"))
            if not hits:
                raise FileNotFoundError(
                    f"panel PDF for sample {sample} not found under {PANEL_DIR}; "
                    "run 05_Fig4e_ExtDataFig9_slide_panels.py first")
            args.pdfs.append(str(hits[0]))

    n_cells = args.rows * args.cols
    if len(args.pdfs) != n_cells:
        parser.error(f"exactly {n_cells} input PDFs are required ({args.rows} x {args.cols}), got {len(args.pdfs)}")

    readers: List[PdfReader] = []
    pages: List[PageObject] = []
    for path in args.pdfs:
        reader = PdfReader(path)
        if len(reader.pages) != 1:
            raise ValueError(f"{path}: expected a single-page PDF, found {len(reader.pages)} pages")
        readers.append(reader)
        pages.append(reader.pages[0])

    widths = [float(p.mediabox.width) for p in pages]
    heights = [float(p.mediabox.height) for p in pages]
    cell_w = max(widths)
    cell_h = max(heights)
    gap = args.gap * PT_PER_INCH
    margin = args.margin * PT_PER_INCH

    page_w = 2 * margin + args.cols * cell_w + (args.cols - 1) * gap
    page_h = 2 * margin + args.rows * cell_h + (args.rows - 1) * gap
    out_page = PageObject.create_blank_page(width=page_w, height=page_h)

    for k, page in enumerate(pages):
        r, c = divmod(k, args.cols)
        w, h = widths[k], heights[k]
        cell_x0 = margin + c * (cell_w + gap)              # left edge of the cell
        cell_y1 = page_h - margin - r * (cell_h + gap)     # top edge of the cell (PDF y grows upward)
        x = cell_x0 + (cell_w - w) / 2.0
        if args.align == "top":
            y = cell_y1 - h
        else:
            y = cell_y1 - cell_h + (cell_h - h) / 2.0
        # Shift by the source page's own media-box origin so its lower-left corner lands at (x, y).
        dx = x - float(page.mediabox.left)
        dy = y - float(page.mediabox.bottom)
        out_page.merge_transformed_page(page, Transformation().translate(dx, dy))
        print(f"[OK]   cell ({r},{c}) {Path(args.pdfs[k]).name}  {w / PT_PER_INCH:.1f} x {h / PT_PER_INCH:.1f} in",
              flush=True)

    out_pdf = Path(args.out_pdf)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    writer = PdfWriter()
    writer.add_page(out_page)
    with open(out_pdf, "wb") as fh:
        writer.write(fh)
    print(f"Saved: {out_pdf} (page {page_w / PT_PER_INCH:.1f} x {page_h / PT_PER_INCH:.1f} in)", flush=True)


if __name__ == "__main__":
    main()
