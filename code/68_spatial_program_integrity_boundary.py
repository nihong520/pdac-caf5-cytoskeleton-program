"""Spatial program integrity and PLOD2–VCL boundary analysis for GSE282302.

This is computational support for program co-localization, NOT wet-lab mechanism
validation. GSE282302 has no patient OS in GEO, so complete-program vs survival
is not performed.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

ROOT = Path.cwd()
OUT = ROOT / "results" / "step24_spatial_program_integrity"
OUT.mkdir(parents=True, exist_ok=True)

GENES5 = ["PLOD2", "MSN", "TPM4", "PPFIBP1", "SEC31A"]
OTHER4 = ["MSN", "TPM4", "PPFIBP1", "SEC31A"]
BOUNDARY = ["PLOD2", "VCL"]


def load_spatial_module():
    path = ROOT / "code" / "63_gse282302_spatial_validation.py"
    spec = importlib.util.spec_from_file_location("gse282302", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def patient_id(sample: str) -> str:
    import re

    m = re.search(r"_(C\d+_D\d+)", sample)
    if not m:
        raise ValueError(sample)
    return m.group(1)


def zscore(values: np.ndarray) -> np.ndarray:
    sd = np.std(values, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(values, dtype=float)
    return (values - np.mean(values)) / sd


def gene_expr(normalized, rows_map, gene: str) -> np.ndarray:
    idxs = rows_map.get(gene, [])
    if not idxs:
        return np.zeros(normalized.shape[1], dtype=float)
    return np.asarray(normalized[idxs[0], :].todense()).ravel()


def process_section(mod, sample: dict) -> dict | None:
    from scipy import sparse

    matrix, symbols, barcodes = mod.read_matrix(sample["h5"])
    positions = mod.read_positions(sample["positions"]).reindex(barcodes)
    keep_pos = positions["in_tissue"].fillna(0).astype(int).eq(1).to_numpy()
    total = np.asarray(matrix.sum(axis=0)).ravel()
    detected = np.asarray((matrix > 0).sum(axis=0)).ravel()
    rows_map = mod.gene_rows(symbols)
    mt_rows = [i for s, idxs in rows_map.items() if s.startswith("MT-") for i in idxs]
    if mt_rows:
        mt = np.asarray(matrix[mt_rows, :].sum(axis=0)).ravel()
        pct_mt = np.divide(mt * 100, total, out=np.zeros_like(mt), where=total > 0)
    else:
        pct_mt = np.zeros_like(total)
    keep = keep_pos & (total >= 500) & (detected >= 200) & (pct_mt < 25)
    if keep.sum() < 200:
        return None

    matrix = matrix[:, keep]
    scale = sparse.diags(1e4 / np.maximum(total[keep], 1.0))
    normalized = (matrix @ scale).astype(np.float32)
    normalized.data = np.log1p(normalized.data)

    expr = {g: gene_expr(normalized, rows_map, g) for g in GENES5 + ["VCL"]}
    z = {g: zscore(expr[g]) for g in expr}
    detected_flag = {g: expr[g] > 0 for g in GENES5}
    n_detected = np.sum(np.vstack([detected_flag[g] for g in GENES5]), axis=0)
    other4_mean_z = zscore(np.mean(np.vstack([z[g] for g in OTHER4]), axis=0))
    complete = n_detected == 5
    incomplete = n_detected <= 2

    # Pairwise Spearman among the five genes within the section
    pair_rows = []
    for i, a in enumerate(GENES5):
        for b in GENES5[i + 1 :]:
            rho, p = stats.spearmanr(expr[a], expr[b])
            pair_rows.append(
                {
                    "Sample": sample["sample"],
                    "Patient": patient_id(sample["sample"]),
                    "Gene_A": a,
                    "Gene_B": b,
                    "Spearman_rho": float(rho),
                    "P": float(p),
                    "Spots": int(keep.sum()),
                }
            )

    # PLOD2-high spots (top quartile of detected PLOD2 or overall PLOD2 z)
    pl_high = z["PLOD2"] >= np.quantile(z["PLOD2"], 0.75)
    pl_low = z["PLOD2"] <= np.quantile(z["PLOD2"], 0.25)
    other4_in_high = float(np.mean(other4_mean_z[pl_high])) if pl_high.any() else np.nan
    other4_in_low = float(np.mean(other4_mean_z[pl_low])) if pl_low.any() else np.nan

    # Completeness co-enrichment: mean CAF-like score proxy = mean of 5 z
    program_z = zscore(np.mean(np.vstack([z[g] for g in GENES5]), axis=0))
    complete_mean = float(np.mean(program_z[complete])) if complete.any() else np.nan
    incomplete_mean = float(np.mean(program_z[incomplete])) if incomplete.any() else np.nan

    # PLOD2 vs VCL
    rho_pv, p_pv = stats.spearmanr(expr["PLOD2"], expr["VCL"])
    both_high = (z["PLOD2"] >= 0.5) & (z["VCL"] >= 0.5)
    plod_only = (z["PLOD2"] >= 0.5) & (z["VCL"] < 0)
    vcl_only = (z["VCL"] >= 0.5) & (z["PLOD2"] < 0)

    summary = {
        "Sample": sample["sample"],
        "Patient": patient_id(sample["sample"]),
        "GSM": sample["gsm"],
        "Spots": int(keep.sum()),
        "Mean_pairwise_rho_5gene": float(np.mean([r["Spearman_rho"] for r in pair_rows])),
        "Frac_complete_5gene": float(np.mean(complete)),
        "Frac_incomplete_le2": float(np.mean(incomplete)),
        "Other4_z_in_PLOD2_high": other4_in_high,
        "Other4_z_in_PLOD2_low": other4_in_low,
        "Delta_Other4_PLOD2_high_minus_low": other4_in_high - other4_in_low,
        "Program_z_complete_spots": complete_mean,
        "Program_z_incomplete_spots": incomplete_mean,
        "PLOD2_VCL_Spearman_rho": float(rho_pv),
        "PLOD2_VCL_P": float(p_pv),
        "Frac_both_high": float(np.mean(both_high)),
        "Frac_PLOD2_only": float(np.mean(plod_only)),
        "Frac_VCL_only": float(np.mean(vcl_only)),
        **{f"Detect_pct_{g}": float(np.mean(expr[g] > 0) * 100) for g in GENES5 + ["VCL"]},
    }
    return {"summary": summary, "pairs": pair_rows}


def aggregate_patients(section_df: pd.DataFrame) -> pd.DataFrame:
    num_cols = [
        c
        for c in section_df.columns
        if c not in {"Sample", "Patient", "GSM"} and pd.api.types.is_numeric_dtype(section_df[c])
    ]
    return (
        section_df.groupby("Patient", as_index=False)
        .agg({**{c: "mean" for c in num_cols}, "Sample": "nunique"})
        .rename(columns={"Sample": "Sections"})
    )


def build_figure(patient_df: pd.DataFrame, pair_patient: pd.DataFrame) -> None:
    sns.set_theme(style="ticks", context="paper")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.2))

    # A: mean pairwise rho among 5 genes
    ax = axes[0, 0]
    sns.histplot(patient_df["Mean_pairwise_rho_5gene"], bins=12, color="#2C7FB8", ax=ax)
    ax.axvline(patient_df["Mean_pairwise_rho_5gene"].median(), color="#D55E00", ls="--")
    ax.set_title("A  Patient mean within-spot 5-gene coexpression")
    ax.set_xlabel("Mean pairwise Spearman rho")
    ax.set_ylabel("Patients")

    # B: heatmap of gene-pair rhos
    ax = axes[0, 1]
    mat = (
        pair_patient.pivot_table(
            index="Gene_A", columns="Gene_B", values="Mean_rho", aggfunc="mean"
        )
        .reindex(index=GENES5, columns=GENES5)
    )
    # make symmetric
    for a in GENES5:
        for b in GENES5:
            if pd.isna(mat.loc[a, b]) and a != b:
                if a in mat.index and b in mat.columns:
                    alt = pair_patient[
                        ((pair_patient.Gene_A == b) & (pair_patient.Gene_B == a))
                        | ((pair_patient.Gene_A == a) & (pair_patient.Gene_B == b))
                    ]
                    if len(alt):
                        mat.loc[a, b] = alt["Mean_rho"].mean()
    np.fill_diagonal(mat.values, 1.0)
    sns.heatmap(
        mat,
        annot=True,
        fmt=".2f",
        cmap="RdBu_r",
        center=0,
        vmin=-0.1,
        vmax=0.6,
        ax=ax,
        cbar_kws={"label": "Mean patient rho"},
    )
    ax.set_title("B  Gene–gene coexpression across patients")

    # C: Other4 enrichment in PLOD2-high spots
    ax = axes[0, 2]
    long = patient_df.melt(
        id_vars=["Patient"],
        value_vars=["Other4_z_in_PLOD2_high", "Other4_z_in_PLOD2_low"],
        var_name="Region",
        value_name="Other4_z",
    )
    long["Region"] = long["Region"].map(
        {
            "Other4_z_in_PLOD2_high": "PLOD2-high",
            "Other4_z_in_PLOD2_low": "PLOD2-low",
        }
    )
    sns.boxplot(data=long, x="Region", y="Other4_z", color="#9ECAE1", ax=ax)
    sns.stripplot(data=long, x="Region", y="Other4_z", color="black", size=3, ax=ax)
    delta = patient_df["Delta_Other4_PLOD2_high_minus_low"]
    w = stats.wilcoxon(delta.dropna(), alternative="greater")
    ax.set_title(f"C  Other4 score in PLOD2 regions\nWilcoxon P={w.pvalue:.2e}")
    ax.set_xlabel("")
    ax.set_ylabel("Mean Other4 z-score")

    # D: program completeness fractions
    ax = axes[1, 0]
    sns.boxplot(
        data=patient_df.melt(
            id_vars=["Patient"],
            value_vars=["Frac_complete_5gene", "Frac_incomplete_le2"],
            var_name="Class",
            value_name="Fraction",
        ).assign(
            Class=lambda d: d["Class"].map(
                {
                    "Frac_complete_5gene": "Complete (5/5)",
                    "Frac_incomplete_le2": "Incomplete (≤2/5)",
                }
            )
        ),
        x="Class",
        y="Fraction",
        color="#A1D99B",
        ax=ax,
    )
    ax.set_title("D  Spot completeness of the five-gene program")
    ax.set_xlabel("")
    ax.set_ylabel("Fraction of spots")

    # E: PLOD2–VCL correlation
    ax = axes[1, 1]
    sns.histplot(patient_df["PLOD2_VCL_Spearman_rho"], bins=12, color="#756BB1", ax=ax)
    ax.axvline(0, color="grey", ls=":")
    ax.axvline(patient_df["PLOD2_VCL_Spearman_rho"].median(), color="#D55E00", ls="--")
    pos = int((patient_df["PLOD2_VCL_Spearman_rho"] > 0).sum())
    ax.set_title(
        f"E  PLOD2–VCL spatial coexpression\n{pos}/{len(patient_df)} patients rho>0"
    )
    ax.set_xlabel("Spearman rho")
    ax.set_ylabel("Patients")

    # F: mutual exclusivity fractions
    ax = axes[1, 2]
    frac = patient_df.melt(
        id_vars=["Patient"],
        value_vars=["Frac_both_high", "Frac_PLOD2_only", "Frac_VCL_only"],
        var_name="Niche",
        value_name="Fraction",
    )
    frac["Niche"] = frac["Niche"].map(
        {
            "Frac_both_high": "Both high",
            "Frac_PLOD2_only": "PLOD2-only",
            "Frac_VCL_only": "VCL-only",
        }
    )
    sns.boxplot(data=frac, x="Niche", y="Fraction", color="#FDBB84", ax=ax)
    sns.stripplot(data=frac, x="Niche", y="Fraction", color="black", size=2.8, ax=ax)
    ax.set_title("F  Boundary niches: co-high vs exclusive")
    ax.set_xlabel("")
    ax.set_ylabel("Fraction of spots")
    ax.tick_params(axis="x", rotation=15)

    fig.suptitle(
        "Spatial integrity of the frozen five-gene program and PLOD2–VCL boundary",
        fontweight="bold",
        y=1.01,
    )
    fig.tight_layout()
    out_png = OUT / "spatial_program_integrity_boundary.png"
    out_pdf = OUT / "spatial_program_integrity_boundary.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print("Wrote", out_png)


def main() -> None:
    mod = load_spatial_module()
    samples = mod.list_samples()
    summaries = []
    pairs = []
    for i, sample in enumerate(samples, 1):
        res = process_section(mod, sample)
        if res is None:
            continue
        summaries.append(res["summary"])
        pairs.extend(res["pairs"])
        if i % 20 == 0:
            print(f"Processed {i}/{len(samples)}")

    section_df = pd.DataFrame(summaries)
    pair_df = pd.DataFrame(pairs)
    patient_df = aggregate_patients(section_df)
    pair_patient = (
        pair_df.groupby(["Patient", "Gene_A", "Gene_B"], as_index=False)
        .agg(Mean_rho=("Spearman_rho", "mean"), Sections=("Sample", "nunique"))
    )

    section_df.to_csv(OUT / "section_program_integrity.csv", index=False)
    patient_df.to_csv(OUT / "patient_program_integrity.csv", index=False)
    pair_df.to_csv(OUT / "section_gene_pairwise_rho.csv", index=False)
    pair_patient.to_csv(OUT / "patient_gene_pairwise_rho.csv", index=False)

    # Formal patient-level tests
    delta = patient_df["Delta_Other4_PLOD2_high_minus_low"].dropna()
    w_other = stats.wilcoxon(delta, alternative="greater")
    rho_pv = patient_df["PLOD2_VCL_Spearman_rho"]
    sign_pv = stats.binomtest(int((rho_pv > 0).sum()), len(rho_pv), 0.5, alternative="greater")
    tests = pd.DataFrame(
        [
            {
                "Test": "Other4_z higher in PLOD2-high vs PLOD2-low spots",
                "N_patients": int(delta.shape[0]),
                "Median_effect": float(delta.median()),
                "P": float(w_other.pvalue),
                "Method": "Wilcoxon signed-rank (one-sided)",
            },
            {
                "Test": "PLOD2–VCL Spearman rho > 0 across patients",
                "N_patients": int(len(rho_pv)),
                "Median_effect": float(rho_pv.median()),
                "P": float(sign_pv.pvalue),
                "Method": "Sign test / binomial",
            },
            {
                "Test": "Mean pairwise 5-gene rho > 0 across patients",
                "N_patients": int(len(patient_df)),
                "Median_effect": float(patient_df["Mean_pairwise_rho_5gene"].median()),
                "P": float(
                    stats.binomtest(
                        int((patient_df["Mean_pairwise_rho_5gene"] > 0).sum()),
                        len(patient_df),
                        0.5,
                        alternative="greater",
                    ).pvalue
                ),
                "Method": "Sign test / binomial",
            },
        ]
    )
    tests.to_csv(OUT / "patient_formal_tests.csv", index=False)

    build_figure(patient_df, pair_patient)

    md = OUT / "SPATIAL_PROGRAM_INTEGRITY_CN.md"
    md.write_text(
        "\n".join(
            [
                "# 空间程序完整性与PLOD2–VCL边界分析",
                "",
                "## 边界声明",
                "- 本分析为计算共区域/共表达支持，不是湿实验机制验证。",
                "- GSE282302 GEO记录未提供患者OS，因此不做“完整程序区域 vs 预后”。",
                "- 正式推断以患者为单位聚合切片。",
                "",
                "## 主要结果",
                f"- 患者层面5基因平均两两Spearman中位数 = {patient_df['Mean_pairwise_rho_5gene'].median():.3f}",
                f"- PLOD2高区Other4评分高于低区：中位差 = {delta.median():.3f}, P = {w_other.pvalue:.2e}",
                f"- PLOD2–VCL相关：中位rho = {rho_pv.median():.3f}; {int((rho_pv > 0).sum())}/{len(rho_pv)} 患者为正",
                f"- 完整5/5 spot比例中位数 = {patient_df['Frac_complete_5gene'].median():.3f}",
                "",
                "## 图",
                "- `spatial_program_integrity_boundary.png`",
            ]
        ),
        encoding="utf-8",
    )
    print(tests.to_string(index=False))
    print("Done:", OUT)


if __name__ == "__main__":
    main()
