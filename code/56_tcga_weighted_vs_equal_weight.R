#!/usr/bin/env Rscript
# Sensitivity: equal-weight 5-gene score vs TCGA multivariate-Cox coefficients.
# Coefficients are fit once in TCGA, frozen, then applied to all cohorts.
# Does NOT replace the locked equal-weight primary score.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(GEOquery)
  library(Biobase)
})

root <- normalizePath(getwd())
pdacr_dir <- file.path(root, "data", "raw", "pdacR")
step19 <- file.path(root, "results", "step19_five_gene_locked_validation")
out_dir <- step19
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pkg_tab <- out_dir

genes5 <- c("PLOD2", "MSN", "TPM4", "PPFIBP1", "SEC31A")
caf_genes <- c(
  "FAP", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "PDPN",
  "THY1", "DCN", "LUM", "TAGLN", "PDGFRB", "SPARC"
)

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

zmat <- function(ex, genes) {
  stopifnot(all(genes %in% rownames(ex)))
  z <- t(scale(t(ex[genes, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  z
}

equal_score <- function(ex, genes) {
  as.numeric(scale(colMeans(zmat(ex, genes))))
}

weighted_score <- function(ex, genes, coefs) {
  z <- zmat(ex, genes)
  raw <- as.numeric(coefs[genes] %*% z)
  as.numeric(scale(raw))
}

caf_score_fn <- function(ex) {
  g <- intersect(caf_genes, rownames(ex))
  equal_score(ex, g)
}

prepare_tcga <- function() {
  obj <- readRDS(file.path(pdacr_dir, "TCGA_PAAD.rds"))
  keep <- which(
    as.character(obj$sampInfo$Filter_by) == "whitelist" &
      as.character(obj$sampInfo$Histological.type.by.DCC) ==
        "pancreas-adenocarcinoma ductal type"
  )
  oid <- "TCGA-2J-AABV-01A"
  if (oid %in% colnames(obj$ex)) {
    keep <- setdiff(keep, which(colnames(obj$ex) == oid))
  }
  list(
    ex = collapse_symbols(obj, keep),
    si = as.data.frame(obj$sampInfo[keep, , drop = FALSE]),
    name = "TCGA_PAAD"
  )
}

prepare_moffitt <- function() {
  obj <- readRDS(file.path(pdacr_dir, "Moffitt_GEO_array.rds"))
  keep <- which(tolower(as.character(obj$sampInfo$specimen_type)) == "primary")
  list(
    ex = collapse_symbols(obj, keep),
    si = as.data.frame(obj$sampInfo[keep, , drop = FALSE]),
    name = "Moffitt_GSE71729"
  )
}

prepare_puleo <- function() {
  obj <- readRDS(file.path(pdacr_dir, "Puleo_array.rds"))
  list(
    ex = collapse_symbols(obj, seq_len(ncol(obj$ex))),
    si = as.data.frame(obj$sampInfo),
    name = "Puleo_EMTAB6134"
  )
}

prepare_gse57495 <- function() {
  matrix_file <- file.path(
    root, "data", "raw", "GSE57495", "GSE57495_series_matrix.txt"
  )
  if (!file.exists(matrix_file)) {
    stop("Place GSE57495_series_matrix.txt in data/raw/GSE57495/")
  }
  eset <- getGEO(filename = matrix_file, getGPL = FALSE)
  if (is.list(eset) && !methods::is(eset, "ExpressionSet")) eset <- eset[[1]]
  probe_expr <- exprs(eset)
  pd <- pData(eset)

  gpl_data_file <- file.path(step19, "GPL15048_data.soft")
  gpl_cache <- file.path(step19, "GPL15048_annotation.rds")
  if (file.exists(gpl_cache)) {
    gpl_tab <- readRDS(gpl_cache)
  } else if (file.exists(gpl_data_file)) {
    gpl_lines <- readLines(gpl_data_file, warn = FALSE)
    table_start <- grep("^!platform_table_begin", gpl_lines)
    table_end <- grep("^!platform_table_end", gpl_lines)
    gpl_tab <- read.delim(
      text = gpl_lines[(table_start + 1):(table_end - 1)],
      check.names = FALSE, quote = "", stringsAsFactors = FALSE
    )
    saveRDS(gpl_tab, gpl_cache)
  } else {
    stop("Missing GPL15048 annotation cache")
  }

  id_col <- grep("^ID$|ID_REF|probe", names(gpl_tab), ignore.case = TRUE, value = TRUE)[1]
  sym_col <- grep(
    "^Gene Symbol$|gene.?symbol|^symbol$",
    names(gpl_tab), ignore.case = TRUE, value = TRUE
  )[1]
  clean_symbol <- function(x) {
    x <- toupper(trimws(as.character(x)))
    x[x %in% c("", "---", "NA")] <- NA_character_
    vapply(strsplit(x, "\\s*///\\s*|\\s*//\\s*|;|,"), function(z) {
      z <- trimws(z)
      z <- z[nzchar(z) & z != "---"]
      if (length(z)) z[1] else NA_character_
    }, character(1))
  }
  symbols <- clean_symbol(gpl_tab[[sym_col]])
  names(symbols) <- as.character(gpl_tab[[id_col]])
  mapped <- symbols[rownames(probe_expr)]
  keep <- !is.na(mapped) & nzchar(mapped)
  mapped_expr <- probe_expr[keep, , drop = FALSE]
  mapped_sym <- mapped[keep]
  probe_var <- apply(mapped_expr, 1, var, na.rm = TRUE)
  ord <- order(probe_var, decreasing = TRUE, na.last = TRUE)
  mapped_expr <- mapped_expr[ord, , drop = FALSE]
  mapped_sym <- mapped_sym[ord]
  mapped_expr <- mapped_expr[!duplicated(mapped_sym), , drop = FALSE]
  rownames(mapped_expr) <- mapped_sym[!duplicated(mapped_sym)]
  write.csv(
    mapped_expr[intersect(c(genes5, caf_genes), rownames(mapped_expr)), , drop = FALSE],
    file.path(step19, "GSE57495_five_gene_expr_matrix.csv")
  )

  extract_characteristic <- function(pd, pattern) {
    hits <- vapply(pd, function(x) {
      any(grepl(pattern, as.character(x), ignore.case = TRUE))
    }, logical(1))
    cols <- names(pd)[hits]
    if (!length(cols)) return(rep(NA_character_, nrow(pd)))
    vals <- rep(NA_character_, nrow(pd))
    for (nm in cols) {
      x <- as.character(pd[[nm]])
      hit <- grepl(pattern, x, ignore.case = TRUE)
      vals[hit] <- x[hit]
    }
    vals
  }
  os_raw <- extract_characteristic(pd, "overall survival")
  event_raw <- extract_characteristic(pd, "vital[._ ]?status")
  os_months <- suppressWarnings(as.numeric(sub(".*:\\s*", "", os_raw)))
  event <- ifelse(
    grepl("dead", event_raw, ignore.case = TRUE), 1,
    ifelse(grepl("alive", event_raw, ignore.case = TRUE), 0, NA)
  )
  si <- data.frame(
    survivalA = os_months,
    `censorA.0yes.1no` = event,
    check.names = FALSE
  )
  rownames(si) <- colnames(mapped_expr)
  list(ex = mapped_expr, si = si, name = "GSE57495")
}

os_frame <- function(co, score, caf) {
  d <- data.frame(
    time = suppressWarnings(as.numeric(co$si$survivalA)),
    event = suppressWarnings(
      as.numeric(as.character(co$si[["censorA.0yes.1no"]]))
    ),
    score = score,
    caf = caf
  )
  d[
    is.finite(d$time) & d$time > 0 & d$event %in% c(0, 1) &
      is.finite(d$score) & is.finite(d$caf),
  ]
}

cox_row <- function(d, score_name, cohort, model, formula) {
  fit <- summary(coxph(formula, data = d))
  data.frame(
    Score = score_name,
    Cohort = cohort,
    Model = model,
    N = fit$n,
    Events = fit$nevent,
    HR = unname(fit$conf.int["score", "exp(coef)"]),
    CI_low = unname(fit$conf.int["score", "lower .95"]),
    CI_high = unname(fit$conf.int["score", "upper .95"]),
    P = unname(fit$coefficients["score", "Pr(>|z|)"])
  )
}

meta_fixed <- function(df, score_name, model) {
  beta <- log(df$HR)
  se <- (log(df$CI_high) - log(df$CI_low)) / (2 * 1.96)
  w <- 1 / se^2
  pooled <- sum(w * beta) / sum(w)
  pooled_se <- sqrt(1 / sum(w))
  q <- sum(w * (beta - pooled)^2)
  df_q <- length(beta) - 1
  i2 <- if (df_q > 0 && q > 0) max(0, (q - df_q) / q) * 100 else 0
  data.frame(
    Score = score_name,
    Model = model,
    N_cohorts = nrow(df),
    HR = exp(pooled),
    CI_low = exp(pooled - 1.96 * pooled_se),
    CI_high = exp(pooled + 1.96 * pooled_se),
    P = 2 * pnorm(-abs(pooled / pooled_se)),
    I2 = i2
  )
}

# ---- Fit TCGA coefficients once ----
tcga <- prepare_tcga()
z_tcga <- zmat(tcga$ex, genes5)
d_fit <- data.frame(
  time = suppressWarnings(as.numeric(tcga$si$survivalA)),
  event = suppressWarnings(
    as.numeric(as.character(tcga$si[["censorA.0yes.1no"]]))
  ),
  t(z_tcga)
)
colnames(d_fit)[-(1:2)] <- genes5
d_fit <- d_fit[
  is.finite(d_fit$time) & d_fit$time > 0 & d_fit$event %in% c(0, 1),
]
form <- as.formula(paste("Surv(time, event) ~", paste(genes5, collapse = " + ")))
fit_mv <- coxph(form, data = d_fit)
coefs <- setNames(as.numeric(coef(fit_mv)[genes5]), genes5)

uni <- do.call(rbind, lapply(genes5, function(g) {
  f <- summary(coxph(as.formula(paste("Surv(time, event) ~", g)), data = d_fit))
  data.frame(
    Gene = g,
    Univariable_HR = unname(f$conf.int[g, "exp(coef)"]),
    Univariable_P = unname(f$coefficients[g, "Pr(>|z|)"])
  )
}))
coef_tab <- data.frame(
  Gene = genes5,
  TCGA_multivariable_coef = as.numeric(coefs[genes5]),
  Equal_weight = 0.2,
  stringsAsFactors = FALSE
)
coef_tab <- merge(coef_tab, uni, by = "Gene", sort = FALSE)
write.csv(
  coef_tab,
  file.path(out_dir, "TCGA_weighted_gene_coefficients.csv"),
  row.names = FALSE
)
write.csv(
  coef_tab,
  file.path(pkg_tab, "TCGA_weighted_gene_coefficients.csv"),
  row.names = FALSE
)

cors <- cor(as.matrix(d_fit[, genes5]))
write.csv(
  round(cors, 3),
  file.path(out_dir, "TCGA_five_gene_correlation.csv")
)
write.csv(
  round(cors, 3),
  file.path(pkg_tab, "TCGA_five_gene_correlation.csv")
)
message("Frozen TCGA multivariable coefficients:")
print(coefs)
message("TCGA five-gene correlation:")
print(round(cors, 3))

# ---- Evaluate scores ----
cohorts <- list(
  prepare_tcga(),
  prepare_moffitt(),
  prepare_puleo(),
  prepare_gse57495()
)

rows <- list()
for (co in cohorts) {
  sc_eq <- equal_score(co$ex, genes5)
  sc_wt <- weighted_score(co$ex, genes5, coefs)
  caf <- caf_score_fn(co$ex)
  d_eq <- os_frame(co, sc_eq, caf)
  d_wt <- os_frame(co, sc_wt, caf)
  rows[[length(rows) + 1]] <- cox_row(
    d_eq, "EqualWeight", co$name, "Univariable", Surv(time, event) ~ score
  )
  rows[[length(rows) + 1]] <- cox_row(
    d_eq, "EqualWeight", co$name, "CAF_adjusted",
    Surv(time, event) ~ score + caf
  )
  rows[[length(rows) + 1]] <- cox_row(
    d_wt, "TCGAweighted", co$name, "Univariable", Surv(time, event) ~ score
  )
  rows[[length(rows) + 1]] <- cox_row(
    d_wt, "TCGAweighted", co$name, "CAF_adjusted",
    Surv(time, event) ~ score + caf
  )
}
res <- do.call(rbind, rows)
write.csv(
  res,
  file.path(out_dir, "equal_vs_TCGAweight_by_cohort.csv"),
  row.names = FALSE
)
write.csv(
  res,
  file.path(pkg_tab, "equal_vs_TCGAweight_by_cohort.csv"),
  row.names = FALSE
)

external <- c("Moffitt_GSE71729", "Puleo_EMTAB6134", "GSE57495")
meta <- do.call(rbind, lapply(c("EqualWeight", "TCGAweighted"), function(sc) {
  do.call(rbind, lapply(c("Univariable", "CAF_adjusted"), function(md) {
    d <- res[res$Score == sc & res$Model == md & res$Cohort %in% external, ]
    meta_fixed(d, sc, md)
  }))
}))
write.csv(
  meta,
  file.path(out_dir, "equal_vs_TCGAweight_external_meta.csv"),
  row.names = FALSE
)
write.csv(
  meta,
  file.path(pkg_tab, "equal_vs_TCGAweight_external_meta.csv"),
  row.names = FALSE
)

message("\n=== Per-cohort ===")
print(res)
message("\n=== External meta ===")
print(meta)

sink(file.path(out_dir, "EQUAL_VS_TCGA_WEIGHTED_SUMMARY_CN.md"))
cat("# 等权 vs TCGA多变量Cox加权（敏感性）\n\n")
cat("主评分仍为等权；本分析仅补充。\n\n")
cat("## TCGA冻结系数\n\n")
print(coef_tab)
cat("\n## 外部meta（Moffitt+Puleo+GSE57495）\n\n")
print(meta)
cat("\n## 解读\n\n")
cat("- 主文继续用等权。\n")
cat("- 若加权与等权外部方向一致，说明结论不依赖等权设定。\n")
cat("- 若加权更弱/不稳，也支持等权（相关基因多变量系数易不稳定）。\n")
sink()

message("Done. Tables written to step19 and manuscript package.")
