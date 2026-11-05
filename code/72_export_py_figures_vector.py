#!/usr/bin/env python
"""Rebuild Fig.4 and consolidated S4 as vector PDF+SVG with editable text (Fig.1 rules).

Fig.4 main panels A–F only (Moran's I reported in text).
S4 localization support: UMAP and PLOD2–Other4 coherence.
Incremental C-index is Supplementary Fig. S3 (R export).

Typography: body 8 pt; panel letters bold 10 pt.
Spatial scatter maps may rasterize points for file size; all labels remain text.
"""

from __future__ import annotations

import importlib.util
import re
import shutil
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "manuscript"
MAIN = PKG / "figures" / "main"
SUPP = PKG / "figures" / "supplementary"
SUB_MAIN = PKG / "submission_pack" / "figures" / "main"
SUB_SUPP = PKG / "submission_pack" / "figures" / "supplementary"
OUT = PKG / "figures" / "results_section2_4_ordered"
SC = ROOT / "results" / "step21_GSE155698_singlecell"
SP = ROOT / "results" / "step23_GSE282302_spatial"
INT = ROOT / "results" / "step24_spatial_program_integrity"
for d in (OUT, MAIN, SUPP, SUB_MAIN, SUB_SUPP):
    d.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "font.size": 8,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.titlesize": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
    }
)

CELL_ORDER = [
    "Fibroblast",
    "Pericyte",
    "Endothelial",
    "Epithelial",
    "T/NK",
    "Myeloid",
    "B/Plasma",
    "Acinar",
    "Mast",
]
MODULE_ORDER = ["Fibroblast", "Pericyte", "Endothelial", "Epithelial", "Immune"]
COMPARTMENT_ORDER = [
    "Fibroblast-rich",
    "Vascular-rich",
    "Epithelial-rich",
    "Immune-rich",
    "Other/mixed",
]
COMPARTMENT_COLORS = {
    "Fibroblast-rich": "#3B82C4",
    "Vascular-rich": "#E69F00",
    "Epithelial-rich": "#2CA25F",
    "Immune-rich": "#B2182B",
    "Other/mixed": "#BDBDBD",
}
GENES5 = ["PLOD2", "MSN", "TPM4", "PPFIBP1", "SEC31A"]


def tag(ax, letter: str) -> None:
    ax.text(
        -0.08,
        1.08,
        letter,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        ha="left",
        clip_on=False,
    )


def normalize_svg(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    out = []
    for line in lines:
        if "<text" in line:
            panel = bool(
                re.search(r">[A-G]</text>", line)
                and (
                    "font-weight:bold" in line.replace(" ", "")
                    or "font-weight:700" in line.replace(" ", "")
                    or "font-weight: bold" in line
                    or "font-weight: 700" in line
                )
            )
            size = "10.00pt" if panel else "8.00pt"
            line = re.sub(r"font-size:\s*[0-9.]+px", f"font-size: {size}", line)
            line = re.sub(r"font-size=\s*\"[0-9.]+\"", f'font-size="{size}"', line)
            line = re.sub(r' textLength="[^"]*"', "", line)
            line = re.sub(r" textLength='[^']*'", "", line)
            line = re.sub(r' lengthAdjust="[^"]*"', "", line)
        out.append(line)
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def save_fig(fig, stem: str, also_main: str | None = None, also_supp: str | None = None) -> None:
    pdf = OUT / f"{stem}.pdf"
    svg = OUT / f"{stem}.svg"
    png = OUT / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    normalize_svg(svg)
    if also_main:
        for ext in ("pdf", "svg", "png"):
            shutil.copy2(OUT / f"{stem}.{ext}", MAIN / f"{also_main}.{ext}")
            shutil.copy2(OUT / f"{stem}.{ext}", SUB_MAIN / f"{also_main}.{ext}")
    if also_supp:
        for ext in ("pdf", "svg", "png"):
            shutil.copy2(OUT / f"{stem}.{ext}", SUPP / f"{also_supp}.{ext}")
            shutil.copy2(OUT / f"{stem}.{ext}", SUB_SUPP / f"{also_supp}.{ext}")
    print(f"Wrote {stem}")


def patient_id(series: pd.Series) -> pd.Series:
    return series.str.extract(r"_(C\d+_D\d+)")[0]


def load_fig4_tables():
    correlations = pd.read_csv(SP / "GSE282302_module_correlations.csv")
    compartments = pd.read_csv(SP / "GSE282302_compartment_summary.csv")
    interfaces = pd.read_csv(SP / "GSE282302_interface_summary.csv")
    qc = pd.read_csv(SP / "GSE282302_QC_summary.csv")
    for table in (correlations, compartments, interfaces, qc):
        table["Patient"] = patient_id(table["Sample"])

    patient_cor = (
        correlations.groupby(["Patient", "Module"], as_index=False)
        .agg(Mean_section_rho=("Spearman_rho", "mean"))
    )
    patient_comp = (
        compartments.groupby(["Patient", "Compartment"], as_index=False)
        .agg(Mean_CAF5_z=("Mean_CAF5_z", "mean"))
    )
    patient_interface = (
        interfaces.groupby(["Patient", "Region"], as_index=False)
        .agg(Mean_CAF5_z=("Mean_CAF5_z", "mean"))
    )
    sc_pb = pd.read_csv(SC / "GSE155698_sample_celltype_pseudobulk.csv")
    sc_dot = pd.read_csv(SC / "GSE155698_gene_celltype_dotplot.csv")
    return {
        "patient_cor": patient_cor,
        "patient_comp": patient_comp,
        "patient_interface": patient_interface,
        "section_qc": qc,
        "sc_pb": sc_pb,
        "sc_dot": sc_dot,
    }


def representative_spots():
    spec = importlib.util.spec_from_file_location(
        "gse282302", ROOT / "code" / "63_gse282302_spatial_validation.py"
    )
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    sample = next(
        item
        for item in mod.list_samples()
        if item["sample"] == "GSM8641055_C2_D9_s5"
    )
    result = mod.process_one(sample)
    if result.get("status") != "ok":
        raise RuntimeError(f"Representative section failed: {result}")
    return result["spot_table"]


def build_figure4(tables, spots) -> None:
    sc_pb = tables["sc_pb"]
    sc_dot = tables["sc_dot"].copy()
    patient_cor = tables["patient_cor"]
    patient_comp = tables["patient_comp"]
    patient_interface = tables["patient_interface"]
    # Moran I retained in results text / QC table; not a main Fig.4 panel

    # PDAC tissue only, patient x cell type with >=20 cells
    pb = sc_pb.copy()
    if "Patient" in pb.columns:
        pb = pb[pb["Patient"].astype(str).str.contains("PDAC", case=False, na=False)]
    cell_n = "Cells" if "Cells" in pb.columns else ("n_cells" if "n_cells" in pb.columns else None)
    if cell_n:
        pb = pb[pb[cell_n] >= 20]
    score_col = "Mean_CAF5" if "Mean_CAF5" in pb.columns else "CAF5_score"
    if score_col not in pb.columns:
        for c in pb.columns:
            if "caf5" in c.lower() or "score" in c.lower():
                score_col = c
                break

    fig = plt.figure(figsize=(14.5, 11.2))
    outer = GridSpec(3, 6, figure=fig, hspace=0.40, wspace=0.55)

    ax = fig.add_subplot(outer[0, :3])
    order = [c for c in CELL_ORDER if c in set(pb["CellType"])]
    sns.boxplot(
        data=pb, x="CellType", y=score_col, order=order, color="#9ECAE1", fliersize=0, ax=ax
    )
    sns.stripplot(
        data=pb, x="CellType", y=score_col, order=order, color="grey", size=2.5, alpha=0.65, ax=ax
    )
    ax.axhline(0, color="grey", linewidth=0.8)
    ax.set_title("Single-cell localization by patient pseudobulk", loc="left", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Mean matched-control score")
    ax.tick_params(axis="x", rotation=35)
    tag(ax, "A")

    ax = fig.add_subplot(outer[0, 3:])
    sc_dot["CellType"] = pd.Categorical(sc_dot["CellType"], categories=CELL_ORDER, ordered=True)
    sc_dot["Gene"] = pd.Categorical(sc_dot["Gene"], categories=GENES5, ordered=True)
    sc_dot = sc_dot.sort_values(["CellType", "Gene"])
    scatter = ax.scatter(
        sc_dot["Gene"].cat.codes,
        sc_dot["CellType"].cat.codes,
        s=20 + 320 * sc_dot["Detected_fraction"],
        c=sc_dot["Mean_expression"],
        cmap="Reds",
        edgecolors="grey",
        linewidths=0.25,
    )
    ax.set_xticks(range(len(GENES5)), GENES5, rotation=35, ha="right")
    ax.set_yticks(range(len(CELL_ORDER)), CELL_ORDER)
    ax.invert_yaxis()
    ax.set_title("Member-gene expression by cell type", loc="left", fontweight="bold")
    ax.grid(color="#E5E5E5", linewidth=0.6)
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Mean log-normalized expression")
    tag(ax, "B")

    maps = outer[1, :4].subgridspec(1, 2, wspace=0.12)
    ax1 = fig.add_subplot(maps[0, 0])
    map1 = ax1.scatter(
        spots["Array_col"],
        -spots["Array_row"],
        c=spots["CAF5_z"],
        cmap="RdBu_r",
        vmin=-2,
        vmax=2,
        s=8,
        linewidths=0,
        rasterized=True,
    )
    ax1.set_aspect("equal")
    ax1.set_axis_off()
    ax1.set_title("Representative section: five-gene score", loc="left", fontweight="bold")
    fig.colorbar(map1, ax=ax1, fraction=0.035, pad=0.01)
    tag(ax1, "C")

    ax2 = fig.add_subplot(maps[0, 1])
    map2 = ax2.scatter(
        spots["Array_col"],
        -spots["Array_row"],
        c=spots["Fibroblast_z"],
        cmap="viridis",
        vmin=-2,
        vmax=2,
        s=8,
        linewidths=0,
        rasterized=True,
    )
    ax2.set_aspect("equal")
    ax2.set_axis_off()
    ax2.set_title("Fibroblast module", loc="left", fontweight="bold")
    fig.colorbar(map2, ax=ax2, fraction=0.035, pad=0.01)

    ax = fig.add_subplot(outer[1, 4:])
    sns.boxplot(
        data=patient_cor,
        x="Module",
        y="Mean_section_rho",
        order=MODULE_ORDER,
        color="#9ECAE1",
        fliersize=0,
        ax=ax,
    )
    sns.stripplot(
        data=patient_cor,
        x="Module",
        y="Mean_section_rho",
        order=MODULE_ORDER,
        color="grey",
        size=2.5,
        alpha=0.55,
        ax=ax,
    )
    ax.axhline(0, color="grey", linewidth=0.8)
    ax.set_title("Spatial correlations by patient (n=39)", loc="left", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Mean section Spearman ρ")
    ax.tick_params(axis="x", rotation=30)
    tag(ax, "D")

    ax = fig.add_subplot(outer[2, :3])
    sns.boxplot(
        data=patient_comp,
        x="Compartment",
        y="Mean_CAF5_z",
        order=COMPARTMENT_ORDER,
        hue="Compartment",
        palette=COMPARTMENT_COLORS,
        legend=False,
        fliersize=0,
        ax=ax,
    )
    ax.axhline(0, color="grey", linewidth=0.8)
    ax.set_title("Spatial compartments by patient", loc="left", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Patient mean five-gene z-score")
    ax.tick_params(axis="x", rotation=35)
    tag(ax, "E")

    ax = fig.add_subplot(outer[2, 3:])
    interface = patient_interface.pivot(index="Patient", columns="Region", values="Mean_CAF5_z")[
        ["Fibroblast_interior", "Fibroblast_interface"]
    ].dropna()
    for _, row in interface.iterrows():
        ax.plot(
            [0, 1],
            [row["Fibroblast_interior"], row["Fibroblast_interface"]],
            color="grey",
            alpha=0.22,
            linewidth=0.7,
        )
    ax.boxplot(
        [interface["Fibroblast_interior"], interface["Fibroblast_interface"]],
        positions=[0, 1],
        widths=0.45,
        patch_artist=True,
        boxprops={"facecolor": "#BDD7E7", "edgecolor": "grey"},
        medianprops={"color": "black"},
    )
    ax.set_xticks([0, 1], ["Fibroblast\ninterior", "Fibroblast\nepithelial-adjacent"])
    ax.set_title("Fibroblast boundary analysis (n=39)", loc="left", fontweight="bold")
    ax.set_ylabel("Patient mean five-gene z-score")
    ax.text(0.5, 0.98, "median paired Δ=0.132; P=0.002", transform=ax.transAxes, ha="center", va="top")
    tag(ax, "F")

    fig.suptitle(
        "Cellular and spatial localization of the frozen five-gene program",
        fontsize=10,
        fontweight="bold",
        y=0.995,
    )
    save_fig(fig, "05_Figure4_singlecell_spatial_localization", also_main="Figure4_singlecell_spatial_localization")
    plt.close(fig)


def build_s4() -> None:
    """Localization support: UMAP and GSE282302 program-integrity panels."""
    meta = pd.read_csv(SC / "GSE155698_cell_metadata_scores.csv.gz")
    integ = pd.read_csv(INT / "patient_program_integrity.csv")

    fig = plt.figure(figsize=(11.5, 4.8))
    gs = GridSpec(1, 2, figure=fig, wspace=0.34, bottom=0.22, top=0.88)

    ax = fig.add_subplot(gs[0, 0])
    plot_df = meta
    if len(plot_df) > 40000:
        plot_df = meta.sample(40000, random_state=20260821)
    types = [c for c in CELL_ORDER if c in set(plot_df["CellType"])]
    palette = dict(zip(types, sns.color_palette("tab10", n_colors=len(types))))
    for ct in types:
        sub = plot_df[plot_df["CellType"] == ct]
        ax.scatter(
            sub["UMAP1"],
            sub["UMAP2"],
            s=1.2,
            c=[palette[ct]],
            linewidths=0,
            rasterized=True,
            label=ct,
            alpha=0.7,
        )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("UMAP of quality-controlled tissue cells", loc="left", fontweight="bold")
    ax.legend(markerscale=4, frameon=False, fontsize=6.5, loc="best")
    tag(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    d = integ[["Other4_z_in_PLOD2_low", "Other4_z_in_PLOD2_high"]].dropna()
    for _, row in d.iterrows():
        ax.plot(
            [0, 1],
            [row["Other4_z_in_PLOD2_low"], row["Other4_z_in_PLOD2_high"]],
            color="grey",
            alpha=0.25,
            linewidth=0.7,
        )
    ax.boxplot(
        [d["Other4_z_in_PLOD2_low"], d["Other4_z_in_PLOD2_high"]],
        positions=[0, 1],
        widths=0.45,
        patch_artist=True,
        boxprops={"facecolor": "#BDD7E7", "edgecolor": "grey"},
        medianprops={"color": "black"},
    )
    ax.set_xticks([0, 1], ["PLOD2-low\nOther4", "PLOD2-high\nOther4"])
    ax.set_ylabel("Other4 score")
    ax.set_title("Other4 score in PLOD2 high vs low spots", loc="left", fontweight="bold")
    tag(ax, "B")

    fig.suptitle(
        "Localization support: annotation and program integrity",
        fontsize=10,
        fontweight="bold",
    )
    fig.text(
        0.01,
        0.015,
        textwrap.fill(
            "Fig. S4 | Localization support: annotation and program integrity. "
            "(A) UMAP of quality-controlled GSE155698 tissue cells. "
            "(B) Remaining-four-gene scores in PLOD2-high versus PLOD2-low regions of GSE282302.",
            width=145,
        ),
        fontsize=8,
        va="bottom",
        ha="left",
    )
    save_fig(fig, "06_FigureS4_localization_support", also_supp="FigureS4_localization_support")
    plt.close(fig)

    src = SP / "GSE282302_spatial_validation.pdf"
    if src.exists():
        shutil.copy2(src, OUT / "07a_source_GSE282302_full_diagnostics.pdf")


def write_readme() -> None:
    text = """# Results §2–4 figures (citation order)

True-vector PDF/SVG exports matching Fig.1 rules.

| Order | File stem | Role |
|------:|-----------|------|
| 01 | `01_Figure2_functional_dissection_locking` | Main Fig. 2 |
| 02 | `02_FigureS2_PPI_community_sensitivity` | Supp. Fig. S2 |
| 03 | `03_Figure3_prognostic_validation_robustness` | Main Fig. 3 |
| 04 | `04_FigureS3_model_diagnostics_platform_sensitivity` | Supp. Fig. S3 (incremental C-index) |
| 05 | `05_Figure4_singlecell_spatial_localization` | Main Fig. 4 (A–F) |
| 06 | `06_FigureS4_localization_support` | Supp. Fig. S4 (UMAP, PLOD2–Other4 coherence) |

Former S5/S6 panels that duplicated Fig.4 or were non-essential were retired from the submission set.

Rebuild:
- R: `code/71_export_fig2_fig3_s2_s3_vector.R`
- Python: `code/72_export_py_figures_vector.py`
"""
    (OUT / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    write_readme()
    print("Loading Fig4 tables...")
    tables = load_fig4_tables()
    print("Loading representative Visium section...")
    spots = representative_spots()
    build_figure4(tables, spots)
    build_s4()
    print("PYTHON_VECTOR_EXPORT_COMPLETE")


if __name__ == "__main__":
    main()
