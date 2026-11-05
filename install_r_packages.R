#!/usr/bin/env Rscript

cran_packages <- c(
  "ggplot2",
  "ggrepel",
  "gridExtra",
  "gridGraphics",
  "httr",
  "igraph",
  "jsonlite",
  "metafor",
  "patchwork",
  "survival",
  "svglite",
  "WGCNA"
)

missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioconductor_packages <- c(
  "Biobase",
  "clusterProfiler",
  "GEOquery",
  "org.Hs.eg.db"
)
missing_bioc <- setdiff(bioconductor_packages, rownames(installed.packages()))
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("R dependencies are installed.")
