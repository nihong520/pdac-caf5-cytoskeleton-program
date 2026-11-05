#!/usr/bin/env Rscript
# Statistical strengthening for the frozen five-gene score:
# proportional hazards, nonlinearity, incremental value, random-effects meta,
# and leave-one-cohort-out meta-analysis.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(splines)
})

set.seed(20260821)
root <- normalizePath(getwd())
pdacr_dir <- file.path(root, "data", "raw", "pdacR")
out_dir <- file.path(root, "results", "step20_statistical_strengthening")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

score_set <- function(ex, genes) {
  present <- intersect(genes, rownames(ex))
  stopifnot(length(present) == length(genes))
  z <- t(scale(t(ex[present, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z)))
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
  si <- as.data.frame(obj$sampInfo[keep, , drop = FALSE])
  ex <- collapse_symbols(obj, keep)
  stage_raw <- tolower(as.character(si$AJCC.pathologic.tumor.stage))
  stage <- ifelse(
    grepl("stage iv", stage_raw), 4,
    ifelse(grepl("stage iii", stage_raw), 3,
      ifelse(grepl("iib", stage_raw), 2, 1)
    )
  )
  grade <- suppressWarnings(
    as.numeric(sub("Grade ([1-4]).*", "\\1", si$Grade))
  )
  list(
    ex = ex,
    si = si,
    clinical = data.frame(
      age10 = suppressWarnings(
        as.numeric(as.character(si$Age.at.initial.pathologic.diagnosis))
      ) / 10,
      stage = stage,
      grade = grade
    ),
    name = "TCGA_PAAD",
    terms = c("age10", "stage", "grade")
  )
}

prepare_moffitt <- function() {
  obj <- readRDS(file.path(pdacr_dir, "Moffitt_GEO_array.rds"))
  keep <- which(tolower(as.character(obj$sampInfo$specimen_type)) == "primary")
  si <- as.data.frame(obj$sampInfo[keep, , drop = FALSE])
  ex <- collapse_symbols(obj, keep)
  clean_num <- function(x) {
    x <- as.character(x)
    x[x %in% c("NaN", "", "NA")] <- NA_character_
    suppressWarnings(as.numeric(x))
  }
  margin_raw <- tolower(as.character(si$Margin))
  margin <- ifelse(
    margin_raw == "negative", 0,
    ifelse(margin_raw %in% c("positive", "close"), 1, NA)
  )
  sex_raw <- toupper(as.character(si$sex))
  sex <- ifelse(sex_raw == "F", 0, ifelse(sex_raw == "M", 1, NA))
  list(
    ex = ex,
    si = si,
    clinical = data.frame(
      T = clean_num(si$T),
      N = clean_num(si$N),
      margin = margin,
      sex = sex
    ),
    name = "Moffitt_GSE71729",
    terms = c("T", "N", "margin", "sex")
  )
}

prepare_puleo <- function() {
  obj <- readRDS(file.path(pdacr_dir, "Puleo_array.rds"))
  si <- as.data.frame(obj$sampInfo)
  ex <- collapse_symbols(obj, seq_len(ncol(obj$ex)))
  tnm <- toupper(as.character(si$TNM.tumor.grade))
  t_stage <- suppressWarnings(as.numeric(sub("^T([1-4]).*", "\\1", tnm)))
  n_stage <- suppressWarnings(as.numeric(sub("^T[1-4]N([01]).*", "\\1", tnm)))
  grade <- suppressWarnings(
    as.numeric(sub("^G([1-4]).*", "\\1", si$Tumor.grade))
  )
  margin <- ifelse(
    toupper(as.character(si$Resection.margin)) == "R0", 0,
    ifelse(toupper(as.character(si$Resection.margin)) == "R1", 1, NA)
  )
  sex <- ifelse(
    tolower(as.character(si$Sex)) == "female", 0,
    ifelse(tolower(as.character(si$Sex)) == "male", 1, NA)
  )
  list(
    ex = ex,
    si = si,
    clinical = data.frame(
      T = t_stage,
      N = n_stage,
      grade = grade,
      margin = margin,
      sex = sex
    ),
    name = "Puleo_EMTAB6134",
    terms = c("T", "N", "grade", "margin", "sex")
  )
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
  stopifnot(!anyNA(idx))
  ex <- ex[, idx, drop = FALSE]
  list(
    ex = ex,
    si = data.frame(
      survivalA = d$OS_months,
      `censorA.0yes.1no` = d$Event,
      check.names = FALSE
    ),
    clinical = data.frame(stage = d$Stage),
    name = "GSE57495",
    terms = "stage"
  )
}

make_data <- function(co) {
  time <- suppressWarnings(as.numeric(co$si$survivalA))
  event <- suppressWarnings(
    as.numeric(as.character(co$si[["censorA.0yes.1no"]]))
  )
  d <- cbind(
    data.frame(
      time = time,
      event = event,
      score = score_set(co$ex, genes5),
      caf = score_set(co$ex, intersect(caf_genes, rownames(co$ex)))
    ),
    co$clinical
  )
  needed <- c("time", "event", "score", "caf", co$terms)
  d <- d[
    is.finite(d$time) & d$time > 0 & d$event %in% c(0, 1) &
      complete.cases(d[, needed, drop = FALSE]),
    ,
    drop = FALSE
  ]
  d
}

fml <- function(terms) {
  as.formula(paste("Surv(time, event) ~", paste(terms, collapse = " + ")))
}

c_index <- function(fit) {
  unname(summary(fit)$concordance[1])
}

increment_one <- function(d, terms, family, boot_n = 1000) {
  base_terms <- terms
  if (family == "Clinical_plus_CAF") base_terms <- c(base_terms, "caf")
  f0 <- fml(base_terms)
  f1 <- fml(c(base_terms, "score"))
  fit0 <- coxph(f0, data = d, x = TRUE, y = TRUE)
  fit1 <- coxph(f1, data = d, x = TRUE, y = TRUE)
  ll0 <- as.numeric(logLik(fit0))
  ll1 <- as.numeric(logLik(fit1))
  lrt <- 2 * (ll1 - ll0)
  delta_c <- c_index(fit1) - c_index(fit0)

  boot_delta <- replicate(boot_n, {
    idx <- sample.int(nrow(d), replace = TRUE)
    db <- d[idx, , drop = FALSE]
    tryCatch({
      b0 <- coxph(f0, data = db)
      b1 <- coxph(f1, data = db)
      c_index(b1) - c_index(b0)
    }, error = function(e) NA_real_)
  })
  boot_delta <- boot_delta[is.finite(boot_delta)]
  ci <- quantile(boot_delta, c(0.025, 0.975), na.rm = TRUE)
  data.frame(
    Model_family = family,
    N = nrow(d),
    Events = sum(d$event),
    Base_C = c_index(fit0),
    Base_plus_score_C = c_index(fit1),
    Delta_C = delta_c,
    Delta_C_boot_low = unname(ci[1]),
    Delta_C_boot_high = unname(ci[2]),
    Bootstrap_success = length(boot_delta),
    LRT_chisq = lrt,
    LRT_df = attr(logLik(fit1), "df") - attr(logLik(fit0), "df"),
    LRT_P = pchisq(
      lrt,
      df = attr(logLik(fit1), "df") - attr(logLik(fit0), "df"),
      lower.tail = FALSE
    ),
    AIC_base = AIC(fit0),
    AIC_base_plus_score = AIC(fit1)
  )
}

ph_detail_rows <- list()

diagnostics_one <- function(d, terms, family, cohort) {
  base_terms <- terms
  if (family == "Clinical_plus_CAF") base_terms <- c(base_terms, "caf")
  fit_lin <- coxph(fml(c(base_terms, "score")), data = d, x = TRUE)
  z <- cox.zph(fit_lin, transform = "km")
  ztab <- as.data.frame(z$table)
  ph_detail_rows[[length(ph_detail_rows) + 1]] <<- cbind(
    Cohort = cohort,
    Model_family = family,
    Term = rownames(ztab),
    ztab,
    row.names = NULL
  )
  ph <- data.frame(
    Term = rownames(ztab),
    Rho_or_chisq = ztab[, 1],
    Chisq = if (ncol(ztab) >= 3) ztab[, 2] else ztab[, 1],
    P = ztab[, ncol(ztab)],
    row.names = NULL
  )
  score_p <- ph$P[ph$Term == "score"]
  global_p <- ph$P[ph$Term == "GLOBAL"]

  fit_spline <- coxph(
    fml(c(base_terms, "ns(score, df=3)")),
    data = d
  )
  cmp <- anova(fit_lin, fit_spline, test = "LRT")
  data.frame(
    Model_family = family,
    N = nrow(d),
    Events = sum(d$event),
    PH_score_P = ifelse(length(score_p), score_p, NA_real_),
    PH_global_P = ifelse(length(global_p), global_p, NA_real_),
    Nonlinearity_LRT_chisq = cmp[2, "Chisq"],
    Nonlinearity_df = cmp[2, "Df"],
    Nonlinearity_P = cmp[2, "Pr(>|Chi|)"]
  )
}

cohorts <- list(
  prepare_tcga(), prepare_moffitt(), prepare_puleo(), prepare_gse57495()
)

diag_rows <- list()
inc_rows <- list()
for (co in cohorts) {
  message("Diagnostics: ", co$name)
  d <- make_data(co)
  for (fam in c("Clinical", "Clinical_plus_CAF")) {
    diag_rows[[length(diag_rows) + 1]] <- cbind(
      Cohort = co$name,
      diagnostics_one(d, co$terms, fam, co$name)
    )
    inc_rows[[length(inc_rows) + 1]] <- cbind(
      Cohort = co$name,
      increment_one(d, co$terms, fam, boot_n = 1000)
    )
  }
}
diagnostics <- do.call(rbind, diag_rows)
incremental <- do.call(rbind, inc_rows)
ph_details <- do.call(rbind, ph_detail_rows)

write.csv(
  diagnostics,
  file.path(out_dir, "five_gene_PH_nonlinearity.csv"),
  row.names = FALSE
)
write.csv(
  ph_details,
  file.path(out_dir, "five_gene_PH_term_details.csv"),
  row.names = FALSE
)
write.csv(
  incremental,
  file.path(out_dir, "five_gene_incremental_value.csv"),
  row.names = FALSE
)

# Meta sensitivity based on the already locked clinical models.
mv <- read.csv(file.path(
  root, "results", "step19_five_gene_locked_validation",
  "five_gene_clinical_multivariable.csv"
))
external <- c("Moffitt_GSE71729", "Puleo_EMTAB6134", "GSE57495")
mv <- mv[
  mv$Cohort %in% external &
    mv$Model %in% c("Clinical_adjusted", "Clinical_plus_CAF"),
]

meta_fit <- function(d, method) {
  yi <- log(d$HR)
  sei <- (log(d$CI_high) - log(d$CI_low)) / (2 * 1.96)
  if (method == "Fixed") {
    w <- 1 / sei^2
    est <- sum(w * yi) / sum(w)
    se <- sqrt(1 / sum(w))
    tau2 <- 0
    q <- sum(w * (yi - est)^2)
    i2 <- if (q > 0) max(0, (q - (length(yi) - 1)) / q) * 100 else 0
  } else {
    if (!requireNamespace("metafor", quietly = TRUE)) {
      stop("metafor is required for REML random-effects meta-analysis")
    }
    fit <- metafor::rma.uni(yi = yi, sei = sei, method = "REML")
    est <- as.numeric(fit$b)
    se <- fit$se
    tau2 <- fit$tau2
    i2 <- fit$I2
  }
  data.frame(
    Method = method,
    K = nrow(d),
    HR = exp(est),
    CI_low = exp(est - 1.96 * se),
    CI_high = exp(est + 1.96 * se),
    P = 2 * pnorm(-abs(est / se)),
    Tau2 = tau2,
    I2 = i2
  )
}

meta_rows <- list()
loco_rows <- list()
for (model in unique(mv$Model)) {
  d <- mv[mv$Model == model,]
  for (method in c("Fixed", "Random_REML")) {
    meta_rows[[length(meta_rows) + 1]] <- cbind(
      Model = model,
      meta_fit(d, method)
    )
  }
  for (drop in d$Cohort) {
    dx <- d[d$Cohort != drop,]
    loco_rows[[length(loco_rows) + 1]] <- cbind(
      Model = model,
      Dropped_cohort = drop,
      meta_fit(dx, "Fixed")
    )
  }
}
meta_sensitivity <- do.call(rbind, meta_rows)
loco <- do.call(rbind, loco_rows)

write.csv(
  meta_sensitivity,
  file.path(out_dir, "five_gene_fixed_random_meta.csv"),
  row.names = FALSE
)
write.csv(
  loco,
  file.path(out_dir, "five_gene_leave_one_cohort_out_meta.csv"),
  row.names = FALSE
)

report <- c(
  "# 五基因评分统计学补强",
  "",
  "## 比例风险与非线性",
  paste(capture.output(print(diagnostics, row.names = FALSE)), collapse = "\n"),
  "",
  "## 增量价值",
  paste(capture.output(print(incremental, row.names = FALSE)), collapse = "\n"),
  "",
  "## 固定/随机效应meta",
  paste(
    capture.output(print(meta_sensitivity, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "## 留一队列meta",
  paste(capture.output(print(loco, row.names = FALSE)), collapse = "\n")
)
writeLines(
  report,
  file.path(out_dir, "STATISTICAL_STRENGTHENING_CN.md")
)

print(diagnostics)
print(incremental)
print(meta_sensitivity)
print(loco)
