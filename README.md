# PDAC CAF5 cytoskeleton-remodeling program

Analysis code for the manuscript:

> CAF-anchored co-expression analysis identifies a reproducible five-gene
> cytoskeleton-remodeling program in pancreatic ductal adenocarcinoma

The frozen program contains **PLOD2, MSN, TPM4, PPFIBP1, and SEC31A**. This
repository matches the submitted analysis: TCGA-PAAD, GSE71729, E-MTAB-6134,
GSE57495, GSE183795, GSE155698, and GSE282302.

## Repository contents

- `code/`: download, analysis, validation, and final figure scripts.
- `docs/`: score definition and recorded analysis environment.
- `data/raw/`: local destination for public source data; large files are ignored.
- `results/`: generated intermediate and final analysis outputs; ignored.
- `manuscript/`: generated figure and table copies; ignored.

No raw expression matrix or patient-level source file is committed.

## Environment

- R 4.6.x
- Python 3.11.x
- Python packages: `pip install -r requirements.txt`
- R packages: `Rscript install_r_packages.R`

The versions used for the submitted analysis are recorded in
`docs/R_SESSION_INFO.txt`, `docs/R_PACKAGE_VERSIONS.txt`, and
`requirements.txt`.

Run every command from the repository root.

## Public data

- TCGA-PAAD, GSE71729, and E-MTAB-6134 processed matrices: `pdacR`
  (`code/03_download_pdacr_cohorts.ps1`).
- GSE57495 series matrix and GSE183795 submitter matrix/family SOFT:
  `code/04_download_geo_bulk_inputs.py`.
- GSE155698 tissue counts: Human Cell Atlas Data Portal
  (`code/58_download_gse155698_tissue_hca.py`).
- GSE282302 Visium files: NCBI GEO
  (`code/62_download_gse282302_spatial.py`).

The download scripts write manifests or checksums where applicable.

## Reproduction order

### 1. Download inputs

```powershell
powershell -ExecutionPolicy Bypass -File code/03_download_pdacr_cohorts.ps1
python code/04_download_geo_bulk_inputs.py
python code/58_download_gse155698_tissue_hca.py
python code/62_download_gse282302_spatial.py
```

### 2. CAF-anchored bulk discovery and external validation

```powershell
Rscript code/41_wgcna_figures_cleaned.R
Rscript code/44_yellow317_mcode_functional.R
Rscript code/46_mcode1_four_blocks_prognosis.R
Rscript code/48_prepare_validate_gse57495_five_gene.R
Rscript code/49_five_gene_multivariable_loo.R
Rscript code/50_five_gene_locked_figures.R
Rscript code/54_four_block_community_sensitivity.R
Rscript code/55_walktrap_vs_biology_score_sensitivity.R
Rscript code/56_tcga_weighted_vs_equal_weight.R
Rscript code/57_five_gene_statistical_strengthening.R
Rscript code/69_single_gene_prognosis_proxy.R
Rscript code/74_validate_gse183795_frozen_five_gene.R
```

### 3. Single-cell and spatial localization

```powershell
python code/59_gse155698_five_gene_singlecell.py
python code/63_gse282302_spatial_validation.py
python code/68_spatial_program_integrity_boundary.py
```

### 4. Final figures

```powershell
Rscript code/70_export_fig1_true_vector.R
Rscript code/71_export_fig2_s2_vector.R
Rscript code/75_export_fig3_with_panel_f.R
python code/72_export_py_figures_vector.py
python code/73_build_bmc_additional_files.py
```

The final submitted figure set is Fig. 1–4 and Supplementary Fig. S1–S4.
The final Supplementary Tables S1–S7 are provided in `derived_tables/` and
are assembled into Additional file 1 by `code/73_build_bmc_additional_files.py`.
Analysis outputs use descriptive filenames; the final export scripts above
apply the submitted figure numbering.

## Frozen score

Within each cohort, each of the five genes is z-standardized across samples.
The equal-weight mean is then re-standardized within cohort. The score is not
calculated if a frozen member is missing. See
`docs/FIVE_GENE_SCORE_DEFINITION.md` (English) and
`docs/FIVE_GENE_SCORE_DEFINITION_LOCKED_CN.md` (Chinese).

## Reproducibility notes

- Random seeds are fixed in the relevant scripts.
- Scripts use paths relative to the repository root.
- Run cohort-level analyses only after the corresponding public data download.
- API- or repository-hosted source files may change availability; accession
  identifiers and download URLs are retained in the scripts.
- The code reproduces computational analyses; it does not imply experimental
  validation or causal mechanism.

## Software record

- Project name: PDAC CAF5 cytoskeleton-remodeling program analysis
- Project home page: https://github.com/nihong520/pdac-caf5-cytoskeleton-program
- Archived version: not yet assigned. `CITATION.cff` labels this copy as version 1.0.0, released 2026-11-03. That label is not a DOI. The DOI is added here after this version is archived on Zenodo.
- Operating system: the submitted analysis was run on Windows 11 x64 (build 26200). R and Python scripts use paths relative to the repository root. `code/03_download_pdacr_cohorts.ps1` is a Windows PowerShell script.
- Programming language: R 4.6.0 and Python 3.11.2
- Other requirements: Python packages listed in `requirements.txt`; R packages installed by `install_r_packages.R`. The versions used for the submitted analysis are in `docs/R_PACKAGE_VERSIONS.txt`, `docs/R_SESSION_INFO.txt`, and `requirements.txt`.
- License: MIT. See `LICENSE`.
- Restrictions to use by non-academics: none. The MIT license permits use, copying, modification, and distribution, including for non-academic purposes, provided the copyright and license notices are retained.

## License

MIT License. See `LICENSE`.
