"""Large-cohort spatial validation of the frozen five-gene program in GSE282302."""

from __future__ import annotations

import gzip
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

SEED = 20260821
RNG = np.random.default_rng(SEED)
ROOT = Path.cwd()
RAW = ROOT / "data" / "raw" / "GSE282302_spatial"
OUT = ROOT / "results" / "step23_GSE282302_spatial"
OUT.mkdir(parents=True, exist_ok=True)

GENES5 = ["PLOD2", "MSN", "TPM4", "PPFIBP1", "SEC31A"]
MODULES = {
    "Fibroblast": [
        "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRA", "SPARC"
    ],
    "Endothelial": ["VWF", "KDR", "EMCN", "ENG", "ESAM", "PECAM1"],
    "Pericyte": ["RGS5", "CSPG4", "MCAM", "NOTCH3", "PDGFRB"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1"],
    "Immune": ["PTPRC", "CD3D", "CD79A", "NKG7", "LST1"],
}
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
NEIGHBOR_OFFSETS = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]


def zscore(values: np.ndarray) -> np.ndarray:
    sd = np.std(values, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(values, dtype=float)
    return (values - np.mean(values)) / sd


def list_samples() -> list[dict]:
    samples = []
    for h5_path in sorted(RAW.glob("*_filtered_feature_bc_matrix.h5")):
        stem = h5_path.name.replace("_filtered_feature_bc_matrix.h5", "")
        pos = RAW / f"{stem}_tissue_positions_list.csv.gz"
        if not pos.exists():
            raise FileNotFoundError(pos)
        gsm = stem.split("_", 1)[0]
        samples.append(
            {
                "gsm": gsm,
                "sample": stem,
                "h5": h5_path,
                "positions": pos,
            }
        )
    if len(samples) < 100:
        raise SystemExit(f"Expected ~108 samples, found {len(samples)}")
    return samples


def read_positions(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        table = pd.read_csv(handle, header=None)
    if table.shape[1] == 6:
        table.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row",
            "pxl_col",
        ]
    else:
        raise ValueError(f"Unexpected positions format: {path}")
    table["barcode"] = table["barcode"].astype(str)
    return table.set_index("barcode")


def read_matrix(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as handle:
        root = "matrix" if "matrix" in handle else list(handle.keys())[0]
        group = handle[root]
        data = group["data"][:].astype(np.float32)
        indices = group["indices"][:].astype(np.int32)
        indptr = group["indptr"][:].astype(np.int32)
        shape = tuple(group["shape"][:])
        barcodes = np.array(
            [
                x.decode() if isinstance(x, bytes) else str(x)
                for x in group["barcodes"][:]
            ]
        )
        feature_group = group["features"]
        names = feature_group["name"][:]
        symbols = np.array(
            [x.decode() if isinstance(x, bytes) else str(x) for x in names]
        ).astype(str)
        # CSC: features x barcodes
        from scipy import sparse

        matrix = sparse.csc_matrix((data, indices, indptr), shape=shape)
    return matrix, symbols, barcodes


def gene_rows(symbols: np.ndarray) -> dict[str, list[int]]:
    mapping: dict[str, list[int]] = {}
    for index, symbol in enumerate(symbols):
        mapping.setdefault(symbol.upper(), []).append(index)
    return mapping


def module_score(matrix, rows: dict[str, list[int]], genes: list[str]) -> np.ndarray:
    selected = [index for gene in genes for index in rows.get(gene.upper(), [])]
    if not selected:
        return np.zeros(matrix.shape[1], dtype=float)
    return np.asarray(matrix[selected, :].mean(axis=0)).ravel()


def build_edges(rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    lookup = {
        (int(r), int(c)): i for i, (r, c) in enumerate(zip(rows, cols, strict=True))
    }
    edges = set()
    for index, (r, c) in enumerate(zip(rows, cols, strict=True)):
        for dr, dc in NEIGHBOR_OFFSETS:
            neighbor = lookup.get((int(r + dr), int(c + dc)))
            if neighbor is not None and neighbor != index:
                edges.add(tuple(sorted((index, neighbor))))
    return np.asarray(sorted(edges), dtype=int)


def moran_i(values: np.ndarray, edges: np.ndarray, permutations: int = 199) -> tuple:
    if len(edges) < 10 or np.std(values) == 0:
        return np.nan, np.nan
    centered = values - np.mean(values)
    denominator = np.sum(centered**2)
    n = len(values)
    e = len(edges)
    observed = (n / (2 * e)) * (
        2 * np.sum(centered[edges[:, 0]] * centered[edges[:, 1]]) / denominator
    )
    null = np.empty(permutations)
    for iteration in range(permutations):
        permuted = RNG.permutation(centered)
        null[iteration] = (n / (2 * e)) * (
            2
            * np.sum(permuted[edges[:, 0]] * permuted[edges[:, 1]])
            / denominator
        )
    p = (1 + np.sum(np.abs(null) >= abs(observed))) / (permutations + 1)
    return float(observed), float(p)


def process_one(sample: dict) -> dict:
    matrix, symbols, barcodes = read_matrix(sample["h5"])
    positions = read_positions(sample["positions"])
    positions = positions.reindex(barcodes)
    keep_pos = positions["in_tissue"].fillna(0).astype(int).eq(1).to_numpy()
    total = np.asarray(matrix.sum(axis=0)).ravel()
    detected = np.asarray((matrix > 0).sum(axis=0)).ravel()
    rows_map = gene_rows(symbols)
    mt_rows = [i for s, idxs in rows_map.items() if s.startswith("MT-") for i in idxs]
    if mt_rows:
        mt = np.asarray(matrix[mt_rows, :].sum(axis=0)).ravel()
        pct_mt = np.divide(mt * 100, total, out=np.zeros_like(mt), where=total > 0)
    else:
        pct_mt = np.zeros_like(total)
    keep = keep_pos & (total >= 500) & (detected >= 200) & (pct_mt < 25)
    if keep.sum() < 200:
        return {
            "sample": sample["sample"],
            "gsm": sample["gsm"],
            "status": "low_spots",
            "spots": int(keep.sum()),
        }

    matrix = matrix[:, keep]
    barcodes = barcodes[keep]
    positions = positions.loc[barcodes]
    total = total[keep]
    detected = detected[keep]
    pct_mt = pct_mt[keep]

    # library-size normalize + log1p
    from scipy import sparse

    scale = sparse.diags(1e4 / np.maximum(total, 1.0))
    normalized = (matrix @ scale).astype(np.float32)
    normalized.data = np.log1p(normalized.data)

    gene_z = []
    detection = {}
    for gene in GENES5:
        idxs = rows_map.get(gene, [])
        if not idxs:
            expr = np.zeros(normalized.shape[1])
        else:
            expr = np.asarray(normalized[idxs[0], :].todense()).ravel()
        detection[gene] = float(np.mean(expr > 0) * 100)
        gene_z.append(zscore(expr))
    caf5 = zscore(np.mean(np.vstack(gene_z), axis=0))

    module_z = {}
    for name, genes in MODULES.items():
        module_z[name] = zscore(module_score(normalized, rows_map, genes))

    module_mat = np.column_stack([module_z[name] for name in MODULES])
    winner = np.argmax(module_mat, axis=1)
    maximum = np.max(module_mat, axis=1)
    names = list(MODULES)
    mapping = {
        "Fibroblast": "Fibroblast-rich",
        "Endothelial": "Vascular-rich",
        "Pericyte": "Vascular-rich",
        "Epithelial": "Epithelial-rich",
        "Immune": "Immune-rich",
    }
    compartment = np.array(
        [
            mapping[names[i]] if value > 0 else "Other/mixed"
            for i, value in zip(winner, maximum, strict=True)
        ]
    )

    array_row = positions["array_row"].to_numpy(dtype=int)
    array_col = positions["array_col"].to_numpy(dtype=int)
    edges = build_edges(array_row, array_col)
    moran, moran_p = moran_i(caf5, edges)

    corr_rows = []
    for name in MODULES:
        rho, p = stats.spearmanr(caf5, module_z[name])
        corr_rows.append(
            {
                "Sample": sample["sample"],
                "GSM": sample["gsm"],
                "Module": name,
                "Spearman_rho": rho,
                "Spot_level_P_descriptive": p,
                "Spots": int(keep.sum()),
            }
        )

    compartment_rows = []
    for group in COMPARTMENT_ORDER:
        values = caf5[compartment == group]
        if len(values):
            compartment_rows.append(
                {
                    "Sample": sample["sample"],
                    "GSM": sample["gsm"],
                    "Compartment": group,
                    "Spots": int(len(values)),
                    "Mean_CAF5_z": float(np.mean(values)),
                    "Median_CAF5_z": float(np.median(values)),
                }
            )

    fibro = compartment == "Fibroblast-rich"
    epithelial = compartment == "Epithelial-rich"
    fibro_interface = np.zeros(len(compartment), dtype=bool)
    for first, second in edges:
        if fibro[first] and epithelial[second]:
            fibro_interface[first] = True
        elif epithelial[first] and fibro[second]:
            fibro_interface[second] = True

    interface_rows = []
    for label, mask in [
        ("Fibroblast_interface", fibro_interface),
        ("Fibroblast_interior", fibro & ~fibro_interface),
        ("Epithelial_rich", epithelial),
        ("Fibroblast_rich", fibro),
    ]:
        values = caf5[mask]
        if len(values):
            interface_rows.append(
                {
                    "Sample": sample["sample"],
                    "GSM": sample["gsm"],
                    "Region": label,
                    "Spots": int(len(values)),
                    "Mean_CAF5_z": float(np.mean(values)),
                }
            )

    # keep a compact spot table only for representative plotting later
    spot_table = pd.DataFrame(
        {
            "Sample": sample["sample"],
            "Barcode": barcodes,
            "Array_row": array_row,
            "Array_col": array_col,
            "CAF5_z": caf5,
            "Fibroblast_z": module_z["Fibroblast"],
            "Epithelial_z": module_z["Epithelial"],
            "Compartment": compartment,
        }
    )

    return {
        "sample": sample["sample"],
        "gsm": sample["gsm"],
        "status": "ok",
        "spots": int(keep.sum()),
        "median_umi": float(np.median(total)),
        "median_genes": float(np.median(detected)),
        "median_mt": float(np.median(pct_mt)),
        "moran": moran,
        "moran_p": moran_p,
        "edges": int(len(edges)),
        "detection": detection,
        "correlations": corr_rows,
        "compartments": compartment_rows,
        "interfaces": interface_rows,
        "spot_table": spot_table,
    }


def paired_wilcoxon(table: pd.DataFrame, a: str, b: str, value="Mean_CAF5_z", key="Compartment"):
    pivot = table.pivot_table(index="Sample", columns=key, values=value)
    if a not in pivot.columns or b not in pivot.columns:
        return {"A": a, "B": b, "N": 0}
    pair = pivot[[a, b]].dropna()
    if len(pair) < 5:
        return {"A": a, "B": b, "N": int(len(pair))}
    test = stats.wilcoxon(pair[a], pair[b])
    return {
        "A": a,
        "B": b,
        "N": int(len(pair)),
        "Median_A": float(pair[a].median()),
        "Median_B": float(pair[b].median()),
        "Median_difference": float(np.median(pair[a] - pair[b])),
        "Wilcoxon_P": float(test.pvalue),
        "Direction_A_gt_B": int(np.sum(pair[a] > pair[b])),
    }


def choose_representatives(correlations: pd.DataFrame, results: list[dict]) -> list[str]:
    fibro = correlations[correlations["Module"] == "Fibroblast"].copy()
    fibro = fibro.sort_values("Spearman_rho", ascending=False)
    picks = []
    # high / mid / low fibroblast correlation with enough spots
    for quantile in [0.9, 0.5, 0.2]:
        target = fibro["Spearman_rho"].quantile(quantile)
        candidate = fibro.iloc[(fibro["Spearman_rho"] - target).abs().argsort()[:1]]
        picks.append(candidate["Sample"].iloc[0])
    # ensure unique and present
    picks = list(dict.fromkeys(picks))
    available = {r["sample"]: r for r in results if r.get("status") == "ok"}
    return [p for p in picks if p in available][:3]


def make_figure(result_bundle: dict) -> None:
    correlations = result_bundle["correlations"]
    compartments = result_bundle["compartments"]
    qc = result_bundle["qc"]
    reps = result_bundle["representatives"]
    spot_map = result_bundle["spot_map"]

    sns.set_theme(style="white", context="paper")
    fig = plt.figure(figsize=(14, 11))
    grid = fig.add_gridspec(3, 3, height_ratios=[1.05, 1.05, 1.0])

    for index, sample_name in enumerate(reps):
        spots = spot_map[sample_name]
        ax = fig.add_subplot(grid[index, 0])
        sc = ax.scatter(
            spots["Array_col"],
            -spots["Array_row"],
            c=spots["CAF5_z"],
            cmap="RdBu_r",
            vmin=-2,
            vmax=2,
            s=8,
            linewidths=0,
        )
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.set_title(
            f"{chr(65 + index)}  {sample_name.split('_',1)[1]}: five-gene",
            loc="left",
            weight="bold",
            fontsize=9,
        )
        fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.01)

        ax2 = fig.add_subplot(grid[index, 1])
        sc2 = ax2.scatter(
            spots["Array_col"],
            -spots["Array_row"],
            c=spots["Fibroblast_z"],
            cmap="viridis",
            vmin=-2,
            vmax=2,
            s=8,
            linewidths=0,
        )
        ax2.set_aspect("equal")
        ax2.set_axis_off()
        ax2.set_title("fibroblast module", loc="left", weight="bold", fontsize=9)
        fig.colorbar(sc2, ax=ax2, fraction=0.04, pad=0.01)

    ax = fig.add_subplot(grid[0, 2])
    module_order = ["Fibroblast", "Pericyte", "Endothelial", "Epithelial", "Immune"]
    sns.boxplot(
        data=correlations,
        x="Module",
        y="Spearman_rho",
        order=module_order,
        color="#9ECAE1",
        fliersize=0,
        ax=ax,
    )
    sns.stripplot(
        data=correlations,
        x="Module",
        y="Spearman_rho",
        order=module_order,
        color="grey",
        size=2.2,
        alpha=0.45,
        ax=ax,
    )
    ax.axhline(0, color="grey", linewidth=0.8)
    ax.set_title("D  Section-level correlations (n sections)", loc="left", weight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Spearman ρ with five-gene score")
    ax.tick_params(axis="x", rotation=25)

    ax = fig.add_subplot(grid[1, 2])
    sns.boxplot(
        data=compartments,
        x="Compartment",
        y="Mean_CAF5_z",
        order=COMPARTMENT_ORDER,
        palette=COMPARTMENT_COLORS,
        hue="Compartment",
        legend=False,
        fliersize=0,
        ax=ax,
    )
    ax.axhline(0, color="grey", linewidth=0.8)
    ax.set_title("E  Compartment means across sections", loc="left", weight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Mean five-gene z-score")
    ax.tick_params(axis="x", rotation=30)

    ax = fig.add_subplot(grid[2, 2])
    ax.hist(qc["Moran_I"].dropna(), bins=20, color="#3182BD", edgecolor="white")
    ax.set_title("F  Spatial autocorrelation of five-gene score", loc="left", weight="bold")
    ax.set_xlabel("Moran's I")
    ax.set_ylabel("Sections")

    fig.suptitle(
        "Spatial validation of the frozen five-gene program in GSE282302 Visium PDAC",
        fontsize=13,
        weight="bold",
    )
    fig.tight_layout()
    fig.savefig(OUT / "GSE282302_spatial_validation.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "GSE282302_spatial_validation.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(bundle: dict) -> None:
    correlations = bundle["correlations"]
    comparisons = bundle["comparisons"]
    qc = bundle["qc"]
    detection = bundle["detection"]
    fibro = correlations[correlations["Module"] == "Fibroblast"]["Spearman_rho"]
    lines = [
        "# GSE282302 大样本空间验证",
        "",
        f"- 可用切片数: {len(qc)}",
        f"- 总spots: {int(qc['Spots'].sum())}",
        f"- Fibroblast Spearman 中位数: {fibro.median():.3f} (IQR {fibro.quantile(0.25):.3f}-{fibro.quantile(0.75):.3f})",
        f"- Moran I 中位数: {qc['Moran_I'].median():.3f}; P<0.05 比例: {(qc['Moran_P'] < 0.05).mean():.1%}",
        "",
        "## 切片级配对比较",
        comparisons.to_string(index=False),
        "",
        "## QC摘要",
        qc.describe().to_string(),
        "",
        "## 五基因检出率中位数(%)",
        detection.groupby("Gene")["Percent_detected"].median().to_string(),
    ]
    (OUT / "GSE282302_SPATIAL_VALIDATION_CN.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    print("\n".join(lines[:12]))


def main() -> None:
    samples = list_samples()
    results = []
    for index, sample in enumerate(samples, start=1):
        print(f"[{index}/{len(samples)}] {sample['sample']}", flush=True)
        results.append(process_one(sample))

    ok = [r for r in results if r.get("status") == "ok"]
    if len(ok) < 20:
        raise SystemExit(f"Too few usable sections: {len(ok)}")

    correlations = pd.DataFrame([row for r in ok for row in r["correlations"]])
    compartments = pd.DataFrame([row for r in ok for row in r["compartments"]])
    interfaces = pd.DataFrame([row for r in ok for row in r["interfaces"]])
    detection = pd.DataFrame(
        [
            {
                "Sample": r["sample"],
                "GSM": r["gsm"],
                "Gene": gene,
                "Percent_detected": pct,
            }
            for r in ok
            for gene, pct in r["detection"].items()
        ]
    )
    qc = pd.DataFrame(
        [
            {
                "Sample": r["sample"],
                "GSM": r["gsm"],
                "Spots": r["spots"],
                "Median_UMI": r["median_umi"],
                "Median_genes": r["median_genes"],
                "Median_percent_mt": r["median_mt"],
                "Edges": r["edges"],
                "Moran_I": r["moran"],
                "Moran_P": r["moran_p"],
            }
            for r in ok
        ]
    )
    failed = pd.DataFrame(
        [
            {"Sample": r["sample"], "GSM": r["gsm"], "Status": r["status"], "Spots": r.get("spots")}
            for r in results
            if r.get("status") != "ok"
        ]
    )

    comparisons = pd.DataFrame(
        [
            paired_wilcoxon(compartments, "Fibroblast-rich", "Epithelial-rich"),
            paired_wilcoxon(compartments, "Vascular-rich", "Epithelial-rich"),
            paired_wilcoxon(
                interfaces,
                "Fibroblast_interface",
                "Fibroblast_interior",
                key="Region",
            ),
            paired_wilcoxon(
                interfaces, "Fibroblast_rich", "Epithelial_rich", key="Region"
            ),
        ]
    )

    # module summary
    module_summary = (
        correlations.groupby("Module")["Spearman_rho"]
        .agg(["count", "median", "mean", lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)])
        .rename(columns={"<lambda_0>": "q25", "<lambda_1>": "q75"})
        .reset_index()
    )
    # sign test vs 0 for fibroblast correlation
    fibro = correlations[correlations["Module"] == "Fibroblast"]["Spearman_rho"]
    sign_p = stats.binomtest(int((fibro > 0).sum()), len(fibro), 0.5).pvalue
    module_summary_path_extra = {
        "Fibroblast_positive_sections": int((fibro > 0).sum()),
        "Fibroblast_total_sections": int(len(fibro)),
        "Fibroblast_sign_test_P": float(sign_p),
        "Fibroblast_wilcoxon_vs0_P": float(stats.wilcoxon(fibro).pvalue),
    }

    correlations.to_csv(OUT / "GSE282302_module_correlations.csv", index=False)
    compartments.to_csv(OUT / "GSE282302_compartment_summary.csv", index=False)
    interfaces.to_csv(OUT / "GSE282302_interface_summary.csv", index=False)
    detection.to_csv(OUT / "GSE282302_five_gene_detection.csv", index=False)
    qc.to_csv(OUT / "GSE282302_QC_summary.csv", index=False)
    comparisons.to_csv(OUT / "GSE282302_section_level_comparisons.csv", index=False)
    module_summary.to_csv(OUT / "GSE282302_module_correlation_summary.csv", index=False)
    pd.Series(module_summary_path_extra).to_csv(
        OUT / "GSE282302_fibroblast_sign_tests.csv", header=["value"]
    )
    if len(failed):
        failed.to_csv(OUT / "GSE282302_failed_sections.csv", index=False)

    reps = choose_representatives(correlations, ok)
    spot_map = {r["sample"]: r["spot_table"] for r in ok if r["sample"] in reps}
    # free memory of non-rep spot tables already not stored beyond ok list; drop others
    for r in ok:
        if r["sample"] not in reps:
            r.pop("spot_table", None)

    bundle = {
        "correlations": correlations,
        "compartments": compartments,
        "comparisons": comparisons,
        "qc": qc,
        "detection": detection,
        "representatives": reps,
        "spot_map": spot_map,
    }
    make_figure(bundle)
    write_report(bundle)
    print("GSE282302_SPATIAL_VALIDATION_COMPLETE", flush=True)
    print(module_summary_path_extra)
    print(comparisons.to_string(index=False))


if __name__ == "__main__":
    main()
