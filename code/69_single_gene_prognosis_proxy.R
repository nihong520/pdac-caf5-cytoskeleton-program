#!/usr/bin/env Rscript
# Single-gene Cox across cohorts: proxy for "not PLOD2 alone".
# True cell-type-resolved bi-PRG prognosis is NOT possible in GSE155698 (no OS).

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})

root <- normalizePath(getwd())
genes5 <- c("PLOD2", "MSN", "TPM4", "PPFIBP1", "SEC31A")
pdacr_dir <- file.path(root, "data", "raw", "pdacR")
out_dir <- file.path(root, "results", "step25_single_gene_prognosis_proxy")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

collapse_symbols <- function(obj, keep_samples) {
  ex <- as.matrix(obj$ex[, keep_samples, drop = FALSE])
  symbols <- toupper(trimws(as.character(obj$featInfo$SYMBOL)))
  keep <- !is.na(symbols) & nzchar(symbols)
  ex <- ex[keep, , drop = FALSE]
  symbols <- symbols[keep]
  sums <- rowsum(ex, group = symbols, reorder = FALSE, na.rm = TRUE)
  counts <- as.numeric(table(factor(symbols, levels = rownames(sums))))
  sums / counts
}

prepare_tcga <- function() {
  obj <- readRDS(file.path(pdacr_dir, "TCGA_PAAD.rds"))
  keep <- which(
    as.character(obj$sampInfo$Filter_by) == "whitelist" &
      as.character(obj$sampInfo$Histological.type.by.DCC) ==
        "pancreas-adenocarcinoma ductal type"
  )
  oid <- "TCGA-2J-AABV-01A"
  if (oid %in% colnames(obj$ex)) keep <- setdiff(keep, which(colnames(obj$ex) == oid))
  si <- as.data.frame(obj$sampInfo[keep, , drop = FALSE])
  ex <- collapse_symbols(obj, keep)
  list(ex = ex, si = si, name = "TCGA_PAAD")
}

prepare_moffitt <- function() {
  obj <- readRDS(file.path(pdacr_dir, "Moffitt_GEO_array.rds"))
  keep <- which(tolower(as.character(obj$sampInfo$specimen_type)) == "primary")
  si <- as.data.frame(obj$sampInfo[keep, , drop = FALSE])
  ex <- collapse_symbols(obj, keep)
  list(ex = ex, si = si, name = "Moffitt_GSE71729")
}

prepare_puleo <- function() {
  obj <- readRDS(file.path(pdacr_dir, "Puleo_array.rds"))
  si <- as.data.frame(obj$sampInfo)
  ex <- collapse_symbols(obj, seq_len(ncol(obj$ex)))
  list(ex = ex, si = si, name = "Puleo_EMTAB6134")
}

prepare_gse57495 <- function() {
  expr_path <- file.path(
    root, "results", "step19_five_gene_locked_validation",
    "GSE57495_five_gene_expr_matrix.csv"
  )
  sample_path <- file.path(
    root, "results", "step19_five_gene_locked_validation",
    "GSE57495_five_gene_samples.csv"
  )
  ex <- as.matrix(read.csv(expr_path, row.names = 1, check.names = FALSE))
  d <- read.csv(sample_path, check.names = FALSE)
  idx <- match(d$Sample, colnames(ex))
  ex <- ex[, idx, drop = FALSE]
  list(
    ex = ex,
    si = data.frame(
      survivalA = d$OS_months,
      `censorA.0yes.1no` = d$Event,
      check.names = FALSE
    ),
    name = "GSE57495"
  )
}

cox_one_gene <- function(co, gene) {
  time <- suppressWarnings(as.numeric(co$si$survivalA))
  event <- suppressWarnings(as.numeric(as.character(co$si[["censorA.0yes.1no"]])))
  if (!gene %in% rownames(co$ex)) {
    return(NULL)
  }
  x <- as.numeric(scale(co$ex[gene, ]))
  d <- data.frame(time = time, event = event, x = x)
  d <- d[is.finite(d$time) & d$time > 0 & d$event %in% c(0, 1) & is.finite(d$x), , drop = FALSE]
  fit <- coxph(Surv(time, event) ~ x, data = d)
  s <- summary(fit)
  data.frame(
    Cohort = co$name,
    Gene = gene,
    N = nrow(d),
    Events = sum(d$event),
    HR = unname(s$coefficients[1, "exp(coef)"]),
    CI_low = unname(s$conf.int[1, "lower .95"]),
    CI_high = unname(s$conf.int[1, "upper .95"]),
    P = unname(s$coefficients[1, "Pr(>|z|)"])
  )
}

fixed_meta <- function(df) {
  # Inverse-variance on logHR
  loghr <- log(df$HR)
  se <- (log(df$CI_high) - log(df$CI_low)) / (2 * 1.96)
  w <- 1 / se^2
  m <- sum(w * loghr) / sum(w)
  se_m <- sqrt(1 / sum(w))
  data.frame(
    Gene = df$Gene[1],
    Model = "External_fixed_meta",
    N_cohorts = nrow(df),
    HR = exp(m),
    CI_low = exp(m - 1.96 * se_m),
    CI_high = exp(m + 1.96 * se_m),
    P = 2 * pnorm(abs(m / se_m), lower.tail = FALSE)
  )
}

cohorts <- list(
  prepare_tcga(),
  prepare_moffitt(),
  prepare_puleo(),
  prepare_gse57495()
)

rows <- list()
for (co in cohorts) {
  for (g in genes5) {
    rows[[length(rows) + 1]] <- cox_one_gene(co, g)
  }
}
by_cohort <- do.call(rbind, rows)
external <- by_cohort[by_cohort$Cohort != "TCGA_PAAD", ]
meta <- do.call(rbind, lapply(split(external, external$Gene), fixed_meta))
meta$Gene <- factor(meta$Gene, levels = rev(genes5))

write.csv(by_cohort, file.path(out_dir, "single_gene_cox_by_cohort.csv"), row.names = FALSE)
write.csv(meta, file.path(out_dir, "single_gene_external_meta.csv"), row.names = FALSE)

p <- ggplot(meta, aes(x = HR, y = Gene)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0, color = "#0072B2") +
  geom_point(size = 2.8, color = "#0072B2") +
  scale_x_log10(limits = c(0.85, 1.8)) +
  labs(
    title = "External single-gene associations (Moffitt + Puleo + GSE57495)",
    subtitle = "Proxy for program breadth; not cell-type-resolved bi-PRG",
    x = "Fixed-effect meta HR per 1 SD",
    y = NULL
  ) +
  theme_classic(base_size = 11)

ggsave(
  file.path(out_dir, "single_gene_external_meta.png"),
  p, width = 7.5, height = 4.2, dpi = 300, bg = "white"
)
ggsave(
  file.path(out_dir, "single_gene_external_meta.pdf"),
  p, width = 7.5, height = 4.2, device = cairo_pdf
)

writeLines(
  c(
    "# 单基因预后代理分析（非bi-PRG）",
    "",
    "## 为什么不是真正的bi-PRG",
    "- GSE155698无患者总体生存结局，无法在单细胞细胞类型层面做预后。",
    "- 本分析仅报告各成员在bulk队列中的单基因Cox及外部固定效应meta。",
    "- 结合已有留一基因结果，用于回答“是否完全由PLOD2驱动”。",
    "",
    "## 结果文件",
    "- single_gene_cox_by_cohort.csv",
    "- single_gene_external_meta.csv",
    "- single_gene_external_meta.png"
  ),
  file.path(out_dir, "SINGLE_GENE_PROGNOSIS_PROXY_CN.md")
)

print(meta)
message("Wrote single-gene prognosis proxy outputs")
