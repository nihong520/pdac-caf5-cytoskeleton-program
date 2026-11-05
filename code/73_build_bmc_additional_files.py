#!/usr/bin/env python
"""Build BMC-style Additional files from submission_pack tables/figures."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.utils.dataframe import dataframe_to_rows

ROOT = Path(__file__).resolve().parents[1]
PACK = ROOT / "manuscript" / "submission_pack"
TAB_DIR = ROOT / "derived_tables"
FIG_DIR = PACK / "figures" / "supplementary"
OUT = PACK / "additional_files"
OUT.mkdir(parents=True, exist_ok=True)

# sheet_name (<=31), csv name, short title, description
# Final submitted set: Supplementary Tables S1–S7.
TABLE_SPECS = [
    (
        "Table_S1",
        "SupplementaryTable_S1_CAF12_gene_coverage.csv",
        "CAF 12-gene panel coverage across cohorts",
        "Measurable CAF-panel genes in TCGA, Moffitt, Puleo and GSE57495. Puleo lacks LUM (K_eff=11).",
    ),
    (
        "Table_S2",
        "SupplementaryTable_S2a_community_agreement.csv",
        "Unsupervised community vs biology-block agreement",
        "ARI/NMI between biology-annotated blocks and Louvain/Walktrap/Fast-greedy/Leiden partitions of the 26-gene PPI.",
    ),
    (
        "Table_S3",
        "SupplementaryTable_S3a_biology5_vs_walktrap5_by_cohort.csv",
        "Biology5 vs Walktrap5 Cox by cohort",
        "Cohort-level Cox results for the biology-annotated PLOD2 five-gene score and Walktrap VCL five-gene score.",
    ),
    (
        "Table_S4a",
        "SupplementaryTable_S4a_five_gene_fixed_random_meta.csv",
        "Fixed- and random-effects external meta",
        "Fixed and REML random-effects pooled HRs for clinical and clinical+CAF models.",
    ),
    (
        "Table_S4b",
        "SupplementaryTable_S4b_single_gene_external_meta.csv",
        "External single-gene meta-analysis",
        "Fixed-effect meta-analysis of each frozen member across Moffitt, Puleo and GSE57495.",
    ),
    (
        "Table_S5",
        "SupplementaryTable_S5_celltype_spatial_marker_panels.csv",
        "Single-cell and spatial marker gene lists",
        "Marker genes used for GSE155698 cell-type scores and for GSE282302 Visium compartment scores.",
    ),
    (
        "Table_S6",
        "SupplementaryTable_S6_GSE155698_localization_tests.csv",
        "GSE155698 fibroblast localization tests",
        "Paired Wilcoxon tests of fibroblast versus other cell types (BH-adjusted).",
    ),
    (
        "Table_S7",
        "SupplementaryTable_S7_spatial_patient_inference.csv",
        "Spatial inference in GSE282302",
        "Section-level Moran's I, five patient-level module correlations, and paired Wilcoxon statistics for Figure 4C–F.",
    ),
]

FIG_ORDER = [
    "FigureS1_WGCNA_construction_QC.pdf",
    "FigureS2_PPI_community_sensitivity.pdf",
    "FigureS3_model_diagnostics_platform_sensitivity.pdf",
    "FigureS4_localization_support.pdf",
]

HEADER_FILL = PatternFill("solid", fgColor="1F4E79")
HEADER_FONT = Font(name="Arial", bold=True, color="FFFFFF", size=10)
BODY_FONT = Font(name="Arial", size=9)
TITLE_FONT = Font(name="Arial", bold=True, size=12)
NOTE_FONT = Font(name="Arial", italic=True, size=9)
THIN = Border(
    left=Side(style="thin", color="D9D9D9"),
    right=Side(style="thin", color="D9D9D9"),
    top=Side(style="thin", color="D9D9D9"),
    bottom=Side(style="thin", color="D9D9D9"),
)


def autosize(ws, max_width: int = 36) -> None:
    for col in ws.columns:
        letter = get_column_letter(col[0].column)
        width = 10
        for cell in col:
            if cell.value is None:
                continue
            width = max(width, min(max_width, len(str(cell.value)) + 2))
        ws.column_dimensions[letter].width = width


def write_readme(wb: Workbook) -> None:
    ws = wb.create_sheet("README", 0)
    ws["A1"] = "Additional file 1 | Supplementary Tables"
    ws["A1"].font = Font(name="Arial", bold=True, size=14)
    ws["A2"] = (
        "CAF-anchored five-gene cytoskeleton-remodeling program in PDAC "
        "(BMC Genomics submission pack)."
    )
    ws["A2"].font = NOTE_FONT
    ws["A3"] = (
        "Tables S1–S7. Patient-level raw score dumps and the non-significant CAF-state "
        "comparison were omitted; Table S5 lists prespecified single-cell and Visium marker genes."
    )
    ws["A3"].font = NOTE_FONT

    headers = ["Sheet", "Table ID", "Title", "Description", "Source CSV"]
    for j, h in enumerate(headers, 1):
        cell = ws.cell(5, j, h)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(wrap_text=True, vertical="center")

    for i, (sheet, csv_name, title, desc) in enumerate(TABLE_SPECS, 6):
        table_id = sheet.replace("Table_", "Table ")
        values = [sheet, table_id, title, desc, csv_name]
        for j, v in enumerate(values, 1):
            cell = ws.cell(i, j, v)
            cell.font = BODY_FONT
            cell.alignment = Alignment(wrap_text=True, vertical="top")
            cell.border = THIN

    ws.column_dimensions["A"].width = 12
    ws.column_dimensions["B"].width = 12
    ws.column_dimensions["C"].width = 42
    ws.column_dimensions["D"].width = 70
    ws.column_dimensions["E"].width = 55
    ws.row_dimensions[5].height = 22
    ws.freeze_panes = "A6"


def write_table_sheet(wb: Workbook, sheet: str, csv_name: str, title: str, desc: str) -> None:
    path = TAB_DIR / csv_name
    df = pd.read_csv(path)
    ws = wb.create_sheet(sheet)
    ws["A1"] = sheet.replace("Table_", "Table ")
    ws["A1"].font = TITLE_FONT
    ws["A2"] = title
    ws["A2"].font = Font(name="Arial", bold=True, size=11)
    ws["A3"] = desc
    ws["A3"].font = NOTE_FONT
    ws.merge_cells(start_row=3, start_column=1, end_row=3, end_column=min(8, max(1, df.shape[1])))

    start_row = 5
    for r_idx, row in enumerate(dataframe_to_rows(df, index=False, header=True), start_row):
        for c_idx, value in enumerate(row, 1):
            cell = ws.cell(r_idx, c_idx, value)
            cell.border = THIN
            cell.alignment = Alignment(wrap_text=True, vertical="top")
            if r_idx == start_row:
                cell.fill = HEADER_FILL
                cell.font = HEADER_FONT
            else:
                cell.font = BODY_FONT
    autosize(ws)
    ws.freeze_panes = f"A{start_row + 1}"


def build_tables_xlsx() -> Path:
    wb = Workbook()
    default = wb.active
    wb.remove(default)
    write_readme(wb)
    for spec in TABLE_SPECS:
        write_table_sheet(wb, *spec)
    out = OUT / "Additional_file_1_Supplementary_Tables.xlsx"
    wb.save(out)
    return out


FIG_CAPTIONS = {
    "FigureS1_WGCNA_construction_QC.pdf": (
        "Fig. S1. Sample quality control for WGCNA discovery cohort. "
        "Hierarchical clustering of TCGA-PAAD samples based on expression dissimilarity, "
        "with trait annotations (CAF score, OS time, OS event). "
        "(A) Before outlier removal (n=133). "
        "(B) After removing outlier TCGA-2J-AABV-01A, yielding the final analytical cohort (n=132)."
    ),
    "FigureS2_PPI_community_sensitivity.pdf": (
        "Fig. S2. Topology sensitivity of biology-annotated blocks versus Walktrap communities. "
        "The same 26-gene PPI network colored by biology-annotated blocks (A) or Walktrap communities (B). "
        "Walktrap yielded four communities (ARI=0.740, NMI=0.770); black outlines mark boundary genes (PLOD2, FN1, VCL)."
    ),
    "FigureS3_model_diagnostics_platform_sensitivity.pdf": (
        "Fig. S3. Change in Harrell's C-index after adding the five-gene score to clinical "
        "or clinical plus CAF models. Points are cohort-level \u0394C; lines are 1,000-bootstrap "
        "percentile 95% CIs. Filled points indicate a likelihood-ratio test P<0.05."
    ),
    "FigureS4_localization_support.pdf": (
        "Fig. S4. Localization support: annotation and program integrity. "
        "(A) UMAP of quality-controlled GSE155698 tissue cells. "
        "(B) Remaining-four-gene scores in PLOD2-high versus PLOD2-low regions of GSE282302."
    ),
}


def build_figures_pdf() -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.image import imread

    out = OUT / "Additional_file_2_Supplementary_Figures.pdf"
    with PdfPages(out) as pdf:
        for name in FIG_ORDER:
            png = FIG_DIR / name.replace(".pdf", ".png")
            if not png.exists():
                raise FileNotFoundError(png)
            img = imread(png)
            h, w = img.shape[:2]
            fig_w = 8.27
            fig_h = min(11.69, max(8.5, fig_w * h / w + 0.4))
            fig = plt.figure(figsize=(fig_w, fig_h))
            ax = fig.add_axes([0.03, 0.03, 0.94, 0.94])
            ax.imshow(img)
            ax.axis("off")
            pdf.savefig(fig)
            plt.close(fig)
    return out


def write_manifest(xlsx: Path, pdf: Path) -> Path:
    lines = [
        "# Additional files list (paste into manuscript)",
        "",
        "Add under a section titled **Additional files** (BMC format):",
        "",
        f"- **Additional file 1** ({xlsx.name}): Excel workbook of Supplementary Tables S1–S7 "
        "(machine-readable; one sheet per table, plus a README index).",
        f"- **Additional file 2** ({pdf.name}): Combined PDF of Supplementary Figures S1–S4.",
        "",
        "Suggested in-text citations:",
        "- CAF-panel coverage (Additional file 1: Table S1)",
        "- Community sensitivity (Additional file 1: Table S2)",
        "- Walktrap vs biology five-gene comparison (Additional file 1: Table S3)",
        "- Prognostic meta and single-gene meta (Additional file 1: Tables S4a, S4b)",
        "- Cell-type/compartment marker panels (Additional file 1: Table S5)",
        "- Single-cell localization tests (Additional file 1: Table S6)",
        "- Spatial inference (Additional file 1: Table S7)",
        "- Sample QC dendrograms (Additional file 2: Figure S1)",
        "",
        "## File checklist for upload",
        f"1. `{xlsx.name}`",
        f"2. `{pdf.name}`",
        "",
        "Main figures (Figure 1–4) should still be uploaded separately as individual figure files "
        "(PDF/SVG preferred), not inside Additional file 2.",
    ]
    out = OUT / "Additional_files_list_for_manuscript.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def main() -> None:
    xlsx = build_tables_xlsx()
    pdf = build_figures_pdf()
    manifest = write_manifest(xlsx, pdf)
    print("WROTE", xlsx)
    print("WROTE", pdf)
    print("WROTE", manifest)
    print("XLSX_MB", round(xlsx.stat().st_size / 1e6, 3))
    print("PDF_MB", round(pdf.stat().st_size / 1e6, 3))


if __name__ == "__main__":
    main()
