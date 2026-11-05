"""Localize the frozen five-gene program in GSE155698 tissue scRNA-seq.

The analysis uses processed count matrices from HCA, performs per-sample QC,
major-cell-type annotation from canonical marker modules, and patient-level
pseudobulk comparisons to avoid treating cells as independent replicates.
"""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import igraph as ig
import leidenalg
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse, stats
from scipy.io import mmread
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from umap import UMAP


SEED = 20260821
np.random.seed(SEED)
ROOT = Path.cwd()
RAW = ROOT / "data" / "raw" / "GSE155698_tissue"
OUT = ROOT / "results" / "step21_GSE155698_singlecell"
OUT.mkdir(parents=True, exist_ok=True)

GENES5 = ["PLOD2", "MSN", "TPM4", "PPFIBP1", "SEC31A"]
CELL_MARKERS = {
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRA"],
    "Pericyte": ["RGS5", "CSPG4", "MCAM", "NOTCH3", "RBP1", "PDGFRB"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "ENG", "ESAM"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1"],
    "T/NK": ["CD3D", "CD3E", "TRBC1", "IL7R", "NKG7", "GNLY"],
    "Myeloid": ["LST1", "TYROBP", "FCER1G", "CTSS", "CD68", "LYZ"],
    "B/Plasma": ["CD79A", "MS4A1", "CD37", "CD74", "MZB1", "JCHAIN"],
    "Acinar": ["PRSS1", "REG1A", "AMY2A", "CTRB2", "CPA1"],
    "Endocrine": ["CHGA", "CHGB", "INS", "IAPP", "PCSK1N"],
    "Mast": ["TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2"],
    "Schwann": ["S100B", "SOX10", "PLP1", "S100A1", "MPZ"],
}
CAF_SUBTYPES = {
    "myCAF-like": ["ACTA2", "TAGLN", "MYL9", "TPM1", "TPM2", "CNN1"],
    "iCAF-like": ["CXCL12", "CFD", "C7", "DPT", "IL6", "CXCL14"],
    "apCAF-like": ["CD74", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1"],
}
TYPE_ORDER = [
    "Fibroblast",
    "Pericyte",
    "Endothelial",
    "Epithelial",
    "T/NK",
    "Myeloid",
    "B/Plasma",
    "Acinar",
    "Endocrine",
    "Mast",
    "Schwann",
]


def sample_records() -> list[dict[str, Path | str]]:
    records = []
    for matrix_path in sorted(RAW.glob("*matrix.mtx.gz")):
        stem = matrix_path.name[: -len("matrix.mtx.gz")]
        sample = stem.rstrip("_")
        features = RAW / f"{stem}features.tsv.gz"
        barcodes = RAW / f"{stem}barcodes.tsv.gz"
        if not features.exists() or not barcodes.exists():
            raise FileNotFoundError(f"Incomplete matrix triplet for {sample}")
        records.append(
            {
                "sample": sample,
                "tissue": "Adjacent" if sample.startswith("AdjNorm") else "PDAC",
                "matrix": matrix_path,
                "features": features,
                "barcodes": barcodes,
            }
        )
    if len(records) != 18:
        raise RuntimeError(f"Expected 18 tissue libraries, found {len(records)}")
    return records


def read_features(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", header=None, compression="gzip")
    if frame.shape[1] < 2:
        raise ValueError(f"Feature table lacks symbols: {path}")
    return pd.DataFrame(
        {
            "gene_id": frame.iloc[:, 0].astype(str),
            "symbol": frame.iloc[:, 1].astype(str).str.upper(),
        }
    )


def sparse_row_mean(matrix: sparse.spmatrix, rows: list[int]) -> np.ndarray:
    if not rows:
        return np.zeros(matrix.shape[1], dtype=np.float32)
    return np.asarray(matrix[rows, :].mean(axis=0)).ravel().astype(np.float32)


def process_samples(records: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    qc_rows = []
    cell_rows = []
    global_sum = None
    global_sumsq = None
    global_n = 0
    reference_features = None

    cache = OUT / "normalized_cache"
    cache.mkdir(exist_ok=True)

    for index, rec in enumerate(records, start=1):
        sample = str(rec["sample"])
        print(f"PROCESS {index}/{len(records)} {sample}", flush=True)
        features = read_features(Path(rec["features"]))
        if reference_features is None:
            reference_features = features
        elif not np.array_equal(
            reference_features["gene_id"].values,
            features["gene_id"].values,
        ):
            raise ValueError(f"Feature order mismatch in {sample}")

        with gzip.open(Path(rec["matrix"]), "rb") as handle:
            counts = mmread(handle).tocsc().astype(np.float32)
        barcodes = pd.read_csv(
            rec["barcodes"],
            sep="\t",
            header=None,
            compression="gzip",
        ).iloc[:, 0].astype(str)
        if counts.shape != (len(features), len(barcodes)):
            raise ValueError(f"Dimension mismatch for {sample}: {counts.shape}")

        total = np.asarray(counts.sum(axis=0)).ravel()
        detected = np.asarray((counts > 0).sum(axis=0)).ravel()
        mt_rows = features["symbol"].str.startswith("MT-").values
        mt = np.asarray(counts[mt_rows, :].sum(axis=0)).ravel()
        pct_mt = np.divide(mt * 100, total, out=np.zeros_like(mt), where=total > 0)
        keep = (
            (detected >= 250)
            & (detected <= 6000)
            & (total >= 500)
            & (pct_mt < 20)
        )
        x = counts[:, keep].tocsc()
        kept_total = total[keep]
        x = x @ sparse.diags(1e4 / kept_total)
        x.data = np.log1p(x.data)
        x = x.astype(np.float32)
        sparse.save_npz(cache / f"{sample}.npz", x)

        sums = np.asarray(x.sum(axis=1)).ravel()
        squared = x.copy()
        squared.data **= 2
        sumsq = np.asarray(squared.sum(axis=1)).ravel()
        global_sum = sums if global_sum is None else global_sum + sums
        global_sumsq = sumsq if global_sumsq is None else global_sumsq + sumsq
        global_n += x.shape[1]

        qc_rows.append(
            {
                "Sample": sample,
                "Tissue": rec["tissue"],
                "Cells_before": counts.shape[1],
                "Cells_after": x.shape[1],
                "Median_features_after": float(np.median(detected[keep])),
                "Median_counts_after": float(np.median(total[keep])),
                "Median_pct_mt_after": float(np.median(pct_mt[keep])),
            }
        )
        cell_rows.extend(
            {
                "Cell": f"{sample}::{barcode}",
                "Sample": sample,
                "Patient": re.sub(r"_[AB]$", "", sample),
                "Tissue": rec["tissue"],
                "nFeature": int(n_feature),
                "nCount": float(n_count),
                "percent_mt": float(percent),
                "Local_index": int(local),
            }
            for local, (barcode, n_feature, n_count, percent) in enumerate(
                zip(
                    barcodes[keep],
                    detected[keep],
                    total[keep],
                    pct_mt[keep],
                    strict=True,
                )
            )
        )

    means = global_sum / global_n
    variances = np.maximum(global_sumsq / global_n - means**2, 0)
    gene_stats = reference_features.copy()
    gene_stats["mean_logexpr"] = means
    gene_stats["variance_logexpr"] = variances
    gene_stats.to_csv(OUT / "GSE155698_gene_expression_stats.csv", index=False)
    return (
        pd.DataFrame(qc_rows),
        pd.DataFrame(cell_rows),
    )


def bh_adjust(pvalues: pd.Series) -> np.ndarray:
    values = pvalues.fillna(1.0).to_numpy(dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted_ranked = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted


def run_analysis(qc: pd.DataFrame, cells: pd.DataFrame) -> dict:
    cache = OUT / "normalized_cache"
    stats_table = pd.read_csv(OUT / "GSE155698_gene_expression_stats.csv")
    symbols = stats_table["symbol"].astype(str).str.upper().to_numpy()
    symbol_rows: dict[str, list[int]] = {}
    for index, symbol in enumerate(symbols):
        symbol_rows.setdefault(symbol, []).append(index)

    required = set(GENES5)
    for markers in CELL_MARKERS.values():
        required.update(markers)
    for markers in CAF_SUBTYPES.values():
        required.update(markers)
    missing = sorted(gene for gene in required if gene not in symbol_rows)
    (OUT / "missing_marker_genes.txt").write_text("\n".join(missing))

    allowed = (
        (stats_table["mean_logexpr"] > 0.01)
        & ~stats_table["symbol"].str.startswith(("MT-", "RPL", "RPS"))
    )
    variable_order = (
        stats_table.loc[allowed, "variance_logexpr"]
        .sort_values(ascending=False)
        .index[:3000]
        .tolist()
    )
    required_rows = [
        row
        for gene in required
        for row in symbol_rows.get(gene, [])
    ]
    hvg_rows = np.array(
        sorted(set(variable_order).union(required_rows)),
        dtype=int,
    )

    blocks = []
    module_data: dict[str, list[np.ndarray]] = {
        key: [] for key in CELL_MARKERS
    }
    subtype_data: dict[str, list[np.ndarray]] = {
        key: [] for key in CAF_SUBTYPES
    }
    gene_data: dict[str, list[np.ndarray]] = {key: [] for key in GENES5}
    sample_slices = {}
    cursor = 0

    # Matched-expression controls for the frozen score.
    means = stats_table["mean_logexpr"].to_numpy()
    target_rows = [symbol_rows[g][0] for g in GENES5]
    excluded = set(target_rows)
    control_rows = []
    for row in target_rows:
        distance = np.abs(means - means[row])
        order = np.argsort(distance)
        chosen = [
            idx
            for idx in order
            if idx not in excluded
            and not symbols[idx].startswith(("MT-", "RPL", "RPS"))
        ][:50]
        control_rows.extend(chosen)
    control_rows = sorted(set(control_rows))
    signature_parts = []

    for sample in qc["Sample"]:
        x = sparse.load_npz(cache / f"{sample}.npz").tocsc()
        n_cells = x.shape[1]
        sample_slices[sample] = slice(cursor, cursor + n_cells)
        cursor += n_cells
        blocks.append(x[hvg_rows, :].T.tocsr())

        for cell_type, markers in CELL_MARKERS.items():
            rows = [
                row
                for gene in markers
                for row in symbol_rows.get(gene, [])
            ]
            module_data[cell_type].append(sparse_row_mean(x, rows))
        for subtype, markers in CAF_SUBTYPES.items():
            rows = [
                row
                for gene in markers
                for row in symbol_rows.get(gene, [])
            ]
            subtype_data[subtype].append(sparse_row_mean(x, rows))
        for gene in GENES5:
            gene_data[gene].append(
                sparse_row_mean(x, symbol_rows.get(gene, []))
            )
        target_score = sparse_row_mean(x, target_rows)
        control_score = sparse_row_mean(x, control_rows)
        signature_parts.append(target_score - control_score)

    x_hvg = sparse.vstack(blocks).tocsr()
    signature = np.concatenate(signature_parts)
    for key in module_data:
        module_data[key] = np.concatenate(module_data[key])
    for key in subtype_data:
        subtype_data[key] = np.concatenate(subtype_data[key])
    for key in gene_data:
        gene_data[key] = np.concatenate(gene_data[key])

    print(f"PCA/UMAP on {x_hvg.shape[0]} cells", flush=True)
    pcs = TruncatedSVD(
        n_components=30,
        random_state=SEED,
    ).fit_transform(x_hvg)
    # Remove sample-specific PC centroids while retaining the global centroid.
    pcs_corrected = pcs.copy()
    global_centroid = pcs.mean(axis=0)
    for sample, section in sample_slices.items():
        pcs_corrected[section] -= pcs[section].mean(axis=0)
        pcs_corrected[section] += global_centroid
    pcs_corrected = StandardScaler().fit_transform(pcs_corrected)

    neighbors = NearestNeighbors(
        n_neighbors=31,
        metric="euclidean",
        n_jobs=-1,
    ).fit(pcs_corrected)
    distances, indices = neighbors.kneighbors(pcs_corrected)
    rows = np.repeat(np.arange(len(cells)), 30)
    cols = indices[:, 1:].reshape(-1)
    edge_distances = distances[:, 1:].reshape(-1)
    scale = np.median(edge_distances[edge_distances > 0])
    weights = np.exp(-edge_distances / scale)
    graph = ig.Graph(
        n=len(cells),
        edges=list(zip(rows.tolist(), cols.tolist(), strict=True)),
        directed=False,
    )
    graph.es["weight"] = weights.tolist()
    graph.simplify(
        multiple=True,
        loops=True,
        combine_edges={"weight": "max"},
    )
    partition = leidenalg.find_partition(
        graph,
        leidenalg.RBConfigurationVertexPartition,
        weights=graph.es["weight"],
        resolution_parameter=0.5,
        seed=SEED,
    )
    cluster = np.asarray(partition.membership)
    embedding = UMAP(
        n_neighbors=30,
        min_dist=0.3,
        metric="cosine",
        random_state=SEED,
        n_jobs=1,
    ).fit_transform(pcs_corrected)

    score_frame = pd.DataFrame(module_data)
    score_z = pd.DataFrame(
        StandardScaler().fit_transform(score_frame),
        columns=score_frame.columns,
    )
    cluster_scores = score_z.assign(Cluster=cluster).groupby("Cluster").mean()
    cluster_labels = cluster_scores.idxmax(axis=1).to_dict()
    cell_type = np.array([cluster_labels[value] for value in cluster])

    cells = cells.copy()
    cells["Cluster"] = cluster
    cells["CellType"] = cell_type
    cells["UMAP1"] = embedding[:, 0]
    cells["UMAP2"] = embedding[:, 1]
    cells["CAF5_score"] = signature
    for gene in GENES5:
        cells[gene] = gene_data[gene]
    for subtype in CAF_SUBTYPES:
        cells[subtype] = subtype_data[subtype]
    cells.to_csv(OUT / "GSE155698_cell_metadata_scores.csv.gz", index=False)
    cluster_scores.assign(
        Assigned_type=pd.Series(cluster_labels)
    ).to_csv(OUT / "GSE155698_cluster_marker_scores.csv")

    return {
        "cells": cells,
        "cluster_scores": cluster_scores,
        "missing": missing,
    }


def summarize_and_plot(result: dict, qc: pd.DataFrame) -> dict:
    cells = result["cells"]
    qc.to_csv(OUT / "GSE155698_QC_summary.csv", index=False)
    counts = (
        cells.groupby(
            ["Patient", "Sample", "Tissue", "CellType"],
            observed=True,
        )
        .size()
        .rename("Cells")
        .reset_index()
    )
    counts.to_csv(OUT / "GSE155698_celltype_counts.csv", index=False)

    pdac = cells[cells["Tissue"] == "PDAC"].copy()
    pseudo = (
        pdac.groupby(["Patient", "CellType"], observed=True)
        .agg(
            Cells=("Cell", "size"),
            Mean_CAF5=("CAF5_score", "mean"),
            Median_CAF5=("CAF5_score", "median"),
        )
        .reset_index()
    )
    pseudo = pseudo[pseudo["Cells"] >= 20].copy()
    pseudo.to_csv(
        OUT / "GSE155698_sample_celltype_pseudobulk.csv",
        index=False,
    )

    # Paired, patient-level localization tests (two-sided).
    comparisons = []
    pivot = pseudo.pivot(
        index="Patient",
        columns="CellType",
        values="Mean_CAF5",
    )
    if "Fibroblast" in pivot:
        for comparator in [x for x in TYPE_ORDER if x != "Fibroblast"]:
            if comparator not in pivot:
                continue
            pair = pivot[["Fibroblast", comparator]].dropna()
            if len(pair) < 5:
                continue
            test = stats.wilcoxon(
                pair["Fibroblast"],
                pair[comparator],
                zero_method="wilcox",
            )
            comparisons.append(
                {
                    "Comparison": f"Fibroblast vs {comparator}",
                    "Paired_samples": len(pair),
                    "Median_difference": float(
                        np.median(pair["Fibroblast"] - pair[comparator])
                    ),
                    "Wilcoxon_statistic": float(test.statistic),
                    "P": float(test.pvalue),
                }
            )
    comparison_table = pd.DataFrame(comparisons)
    if not comparison_table.empty:
        comparison_table["P_BH"] = bh_adjust(comparison_table["P"])
    comparison_table.to_csv(
        OUT / "GSE155698_fibroblast_localization_tests.csv",
        index=False,
    )

    # Gene dot-plot values averaged at sample level.
    gene_rows = []
    for (patient, cell_type), group in pdac.groupby(
        ["Patient", "CellType"],
        observed=True,
    ):
        if len(group) < 20:
            continue
        for gene in GENES5:
            gene_rows.append(
                {
                    "Patient": patient,
                    "CellType": cell_type,
                    "Gene": gene,
                    "Mean_expression": group[gene].mean(),
                    "Detected_fraction": (group[gene] > 0).mean(),
                }
            )
    gene_sample = pd.DataFrame(gene_rows)
    gene_sample.to_csv(
        OUT / "GSE155698_gene_celltype_sample_summary.csv",
        index=False,
    )
    dot = (
        gene_sample.groupby(["CellType", "Gene"], observed=True)
        .agg(
            Mean_expression=("Mean_expression", "mean"),
            Detected_fraction=("Detected_fraction", "mean"),
            Patients=("Patient", "nunique"),
        )
        .reset_index()
    )
    dot.to_csv(OUT / "GSE155698_gene_celltype_dotplot.csv", index=False)

    # Fibroblast state analysis, restricted to PDAC.
    fibro = pdac[pdac["CellType"] == "Fibroblast"].copy()
    subtype_columns = list(CAF_SUBTYPES)
    if not fibro.empty:
        subtype_z = StandardScaler().fit_transform(fibro[subtype_columns])
        fibro["CAF_subtype"] = np.array(subtype_columns)[
            np.argmax(subtype_z, axis=1)
        ]
    subtype_pseudo = (
        fibro.groupby(["Patient", "CAF_subtype"], observed=True)
        .agg(
            Cells=("Cell", "size"),
            Mean_CAF5=("CAF5_score", "mean"),
        )
        .reset_index()
    )
    subtype_pseudo = subtype_pseudo[subtype_pseudo["Cells"] >= 20]
    subtype_pseudo.to_csv(
        OUT / "GSE155698_CAF_subtype_pseudobulk.csv",
        index=False,
    )

    subtype_tests = []
    subtype_pivot = subtype_pseudo.pivot(
        index="Patient",
        columns="CAF_subtype",
        values="Mean_CAF5",
    )
    for first, second in [
        ("myCAF-like", "iCAF-like"),
        ("myCAF-like", "apCAF-like"),
        ("iCAF-like", "apCAF-like"),
    ]:
        if first not in subtype_pivot or second not in subtype_pivot:
            continue
        pair = subtype_pivot[[first, second]].dropna()
        if len(pair) < 5:
            continue
        test = stats.wilcoxon(pair[first], pair[second])
        subtype_tests.append(
            {
                "Comparison": f"{first} vs {second}",
                "Paired_samples": len(pair),
                "Median_difference": float(
                    np.median(pair[first] - pair[second])
                ),
                "P": float(test.pvalue),
            }
        )
    subtype_test_table = pd.DataFrame(subtype_tests)
    if not subtype_test_table.empty:
        subtype_test_table["P_BH"] = bh_adjust(subtype_test_table["P"])
    subtype_test_table.to_csv(
        OUT / "GSE155698_CAF_subtype_tests.csv",
        index=False,
    )

    # Figure.
    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.2))
    ax = axes[0, 0]
    plot_cells = cells
    if len(plot_cells) > 35000:
        plot_cells = plot_cells.sample(35000, random_state=SEED)
    palette = dict(
        zip(
            TYPE_ORDER,
            sns.color_palette("tab20", n_colors=len(TYPE_ORDER)),
            strict=True,
        )
    )
    for cell_type in TYPE_ORDER:
        group = plot_cells[plot_cells["CellType"] == cell_type]
        if group.empty:
            continue
        ax.scatter(
            group["UMAP1"],
            group["UMAP2"],
            s=2,
            alpha=0.55,
            color=palette[cell_type],
            label=cell_type,
            rasterized=True,
        )
    ax.set_title("A  GSE155698 tissue cells: major cell types", loc="left", weight="bold")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=7,
        markerscale=3,
    )

    ax = axes[0, 1]
    present_order = [x for x in TYPE_ORDER if x in pseudo["CellType"].unique()]
    sns.boxplot(
        data=pseudo,
        x="CellType",
        y="Mean_CAF5",
        order=present_order,
        color="#d9d9d9",
        fliersize=0,
        ax=ax,
    )
    sns.stripplot(
        data=pseudo,
        x="CellType",
        y="Mean_CAF5",
        order=present_order,
        hue="Patient",
        palette="husl",
        size=3,
        alpha=0.8,
        legend=False,
        ax=ax,
    )
    ax.axhline(0, color="grey", linewidth=0.6)
    ax.set_title(
        "B  Five-gene program by cell type (PDAC patient pseudobulk)",
        loc="left",
        weight="bold",
    )
    ax.set_xlabel("")
    ax.set_ylabel("Mean matched-control score")
    ax.tick_params(axis="x", rotation=45)

    ax = axes[1, 0]
    dot["CellType"] = pd.Categorical(
        dot["CellType"],
        categories=present_order,
        ordered=True,
    )
    dot["Gene"] = pd.Categorical(dot["Gene"], categories=GENES5, ordered=True)
    dot = dot.dropna(subset=["CellType"]).sort_values(["CellType", "Gene"])
    x_map = {gene: i for i, gene in enumerate(GENES5)}
    y_map = {name: i for i, name in enumerate(reversed(present_order))}
    scatter = ax.scatter(
        [x_map[str(x)] for x in dot["Gene"]],
        [y_map[str(y)] for y in dot["CellType"]],
        s=20 + 450 * dot["Detected_fraction"],
        c=dot["Mean_expression"],
        cmap="Reds",
        edgecolors="grey",
        linewidths=0.3,
    )
    ax.set_xticks(range(len(GENES5)), GENES5, rotation=45, ha="right")
    ax.set_yticks(
        range(len(present_order)),
        list(reversed(present_order)),
    )
    ax.set_title(
        "C  Gene expression by cell type (sample-averaged)",
        loc="left",
        weight="bold",
    )
    fig.colorbar(scatter, ax=ax, label="Mean log-normalized expression")

    ax = axes[1, 1]
    subtype_order = ["myCAF-like", "iCAF-like", "apCAF-like"]
    if not subtype_pseudo.empty:
        sns.boxplot(
            data=subtype_pseudo,
            x="CAF_subtype",
            y="Mean_CAF5",
            order=subtype_order,
            color="#f2c9d4",
            fliersize=0,
            ax=ax,
        )
        sns.stripplot(
            data=subtype_pseudo,
            x="CAF_subtype",
            y="Mean_CAF5",
            order=subtype_order,
        hue="Patient",
            palette="husl",
            size=4,
            legend=False,
            ax=ax,
        )
    ax.axhline(0, color="grey", linewidth=0.6)
    ax.set_title(
        "D  Program across CAF-like states (PDAC patient pseudobulk)",
        loc="left",
        weight="bold",
    )
    ax.set_xlabel("")
    ax.set_ylabel("Mean matched-control score")
    ax.tick_params(axis="x", rotation=20)

    fig.suptitle(
        "Single-cell localization of the frozen five-gene program",
        fontsize=14,
        weight="bold",
    )
    fig.tight_layout()
    png = OUT / "GSE155698_singlecell_localization.png"
    pdf = OUT / "GSE155698_singlecell_localization.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    return {
        "qc": qc,
        "counts": counts,
        "pseudo": pseudo,
        "comparisons": comparison_table,
        "subtype_pseudo": subtype_pseudo,
        "subtype_tests": subtype_test_table,
    }


def write_report(summary: dict, missing: list[str]) -> None:
    lines = [
        "# GSE155698单细胞定位：冻结五基因程序",
        "",
        "## 数据与原则",
        "- 仅使用组织样本（PDAC及邻近组织），不使用PBMC。",
        "- 细胞层面仅用于可视化和注释；正式比较以样本×细胞类型pseudobulk为单位。",
        "- 评分为五基因平均log-normalized表达减去表达匹配对照基因均值。",
        "",
        "## QC",
        summary["qc"].to_string(index=False),
        "",
        "## PDAC样本层面成纤维细胞定位检验",
        summary["comparisons"].to_string(index=False),
        "",
        "## CAF-like状态比较",
        summary["subtype_tests"].to_string(index=False),
        "",
        "## 缺失标志基因",
        ", ".join(missing) if missing else "None",
        "",
        "## 边界",
        "- 主要细胞类型与CAF-like状态为标志基因模块注释，不等同于组织学或谱系追踪。",
        "- 不将单个细胞视为独立生物学重复。",
    ]
    (OUT / "GSE155698_SINGLECELL_LOCALIZATION_CN.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def main() -> None:
    records = sample_records()
    qc, cells = process_samples(records)
    result = run_analysis(qc, cells)
    summary = summarize_and_plot(result, qc)
    write_report(summary, result["missing"])
    print(summary["comparisons"].to_string(index=False))
    print(summary["subtype_tests"].to_string(index=False))
    print("SINGLECELL_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
