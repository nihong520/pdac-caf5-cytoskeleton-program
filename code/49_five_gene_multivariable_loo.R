#!/usr/bin/env Rscript
# Locked 5-gene validation: clinical multivariable models and leave-one-gene-out.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(survival))

project_root <- normalizePath(getwd())
pdacr_dir <- file.path(project_root, "data", "raw", "pdacR")
out_dir <- file.path(project_root, "results", "step19_five_gene_locked_validation")
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
  rownames(si) <- colnames(ex)
  stage_raw <- tolower(as.character(si$AJCC.pathologic.tumor.stage))
  stage <- ifelse(
    grepl("stage iv", stage_raw), 4,
    ifelse(
      grepl("stage iii", stage_raw), 3,
      ifelse(grepl("iib", stage_raw), 2, 1)
    )
  )
  grade <- suppressWarnings(as.numeric(sub("Grade ([1-4]).*", "\\1", si$Grade)))
  clinical <- data.frame(
    age10 = suppressWarnings(
      as.numeric(as.character(si$Age.at.initial.pathologic.diagnosis))
    ) / 10,
    stage = stage,
    grade = grade,
    subtype = factor(as.character(si$Moffitt.Clusters))
  )
  list(
    ex = ex, si = si, clinical = clinical, name = "TCGA_PAAD",
    clinical_terms = c("age10", "stage", "grade")
  )
}

prepare_moffitt <- function() {
  obj <- readRDS(file.path(pdacr_dir, "Moffitt_GEO_array.rds"))
  keep <- which(tolower(as.character(obj$sampInfo$specimen_type)) == "primary")
  si <- as.data.frame(obj$sampInfo[keep, , drop = FALSE])
  ex <- collapse_symbols(obj, keep)
  rownames(si) <- colnames(ex)
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
  clinical <- data.frame(
    T = clean_num(si$T),
    N = clean_num(si$N),
    margin = margin,
    sex = sex,
    subtype = factor(as.character(si$MoffittTumor))
  )
  list(
    ex = ex, si = si, clinical = clinical, name = "Moffitt_GSE71729",
    clinical_terms = c("T", "N", "margin", "sex")
  )
}

prepare_puleo <- function() {
  obj <- readRDS(file.path(pdacr_dir, "Puleo_array.rds"))
  si <- as.data.frame(obj$sampInfo)
  ex <- collapse_symbols(obj, seq_len(ncol(obj$ex)))
  rownames(si) <- colnames(ex)
  tnm <- toupper(as.character(si$TNM.tumor.grade))
  t_stage <- suppressWarnings(as.numeric(sub("^T([1-4]).*", "\\1", tnm)))
  n_stage <- suppressWarnings(as.numeric(sub("^T[1-4]N([01]).*", "\\1", tnm)))
  grade <- suppressWarnings(as.numeric(sub("^G([1-4]).*", "\\1", si$Tumor.grade)))
  margin <- ifelse(
    toupper(as.character(si$Resection.margin)) == "R0", 0,
    ifelse(toupper(as.character(si$Resection.margin)) == "R1", 1, NA)
  )
  sex <- ifelse(
    tolower(as.character(si$Sex)) == "female", 0,
    ifelse(tolower(as.character(si$Sex)) == "male", 1, NA)
  )
  clinical <- data.frame(
    T = t_stage,
    N = n_stage,
    grade = grade,
    margin = margin,
    sex = sex,
    subtype = factor(as.character(si$WholeTumorClass))
  )
  list(
    ex = ex, si = si, clinical = clinical, name = "Puleo_EMTAB6134",
    clinical_terms = c("T", "N", "grade", "margin", "sex")
  )
}

prepare_gse57495 <- function() {
  ex <- readRDS(file.path(out_dir, "GSE57495_expr_genes.rds"))
  d <- read.csv(
    file.path(out_dir, "GSE57495_five_gene_samples.csv"),
    stringsAsFactors = FALSE
  )
  idx <- match(d$Sample, colnames(ex))
  ex <- ex[, idx, drop = FALSE]
  si <- data.frame(
    survivalA = d$OS_months,
    censorA.0yes.1no = d$Event
  )
  clinical <- data.frame(stage = d$Stage, subtype = factor(NA))
  list(
    ex = ex, si = si, clinical = clinical, name = "GSE57495",
    clinical_terms = "stage"
  )
}

cohorts <- list(
  prepare_tcga(), prepare_moffitt(), prepare_puleo(), prepare_gse57495()
)

extract_score <- function(fit, model, cohort, complete_n) {
  sm <- summary(fit)
  data.frame(
    Cohort = cohort,
    Model = model,
    N = complete_n,
    Events = sm$nevent,
    HR = unname(sm$conf.int["score", "exp(coef)"]),
    CI_low = unname(sm$conf.int["score", "lower .95"]),
    CI_high = unname(sm$conf.int["score", "upper .95"]),
    P = unname(sm$coefficients["score", "Pr(>|z|)"])
  )
}

multivariable_rows <- list()
for (co in cohorts) {
  time <- suppressWarnings(as.numeric(co$si$survivalA))
  event <- suppressWarnings(
    as.numeric(as.character(co$si[["censorA.0yes.1no"]]))
  )
  d <- cbind(
    data.frame(
      time = time,
      event = event,
      score = score_set(co$ex, genes5),
      caf = score_set(co$ex, caf_genes)
    ),
    co$clinical
  )
  base_ok <- is.finite(d$time) & d$time > 0 & d$event %in% c(0, 1)

  clinical_formula <- as.formula(paste(
    "Surv(time, event) ~ score +",
    paste(co$clinical_terms, collapse = " + ")
  ))
  needed <- c("time", "event", "score", co$clinical_terms)
  dc <- d[base_ok & complete.cases(d[, needed, drop = FALSE]), ]
  fit_clinical <- coxph(clinical_formula, data = dc)
  multivariable_rows[[length(multivariable_rows) + 1]] <- extract_score(
    fit_clinical, "Clinical_adjusted", co$name, nrow(dc)
  )

  clinical_caf_formula <- update(clinical_formula, . ~ . + caf)
  needed_caf <- c(needed, "caf")
  dcc <- d[base_ok & complete.cases(d[, needed_caf, drop = FALSE]), ]
  fit_clinical_caf <- coxph(clinical_caf_formula, data = dcc)
  multivariable_rows[[length(multivariable_rows) + 1]] <- extract_score(
    fit_clinical_caf, "Clinical_plus_CAF", co$name, nrow(dcc)
  )

  # Expression-derived subtype sensitivity analysis where available.
  if ("subtype" %in% names(d) && nlevels(droplevels(d$subtype)) >= 2) {
    ds <- d[base_ok & complete.cases(d[, c("score", "subtype")]), ]
    if (nrow(ds) >= 40) {
      fit_subtype <- coxph(Surv(time, event) ~ score + subtype, data = ds)
      multivariable_rows[[length(multivariable_rows) + 1]] <- extract_score(
        fit_subtype, "Subtype_adjusted", co$name, nrow(ds)
      )
    }
  }
}
multivariable <- do.call(rbind, multivariable_rows)
write.csv(
  multivariable,
  file.path(out_dir, "five_gene_clinical_multivariable.csv"),
  row.names = FALSE
)

# Leave-one-gene-out: full score plus five four-gene variants. No variant selection.
variants <- c(
  list(Full5 = genes5),
  setNames(lapply(genes5, function(g) setdiff(genes5, g)), paste0("Drop_", genes5))
)
loo_rows <- list()
for (variant in names(variants)) {
  for (co in cohorts) {
    time <- suppressWarnings(as.numeric(co$si$survivalA))
    event <- suppressWarnings(
      as.numeric(as.character(co$si[["censorA.0yes.1no"]]))
    )
    score <- score_set(co$ex, variants[[variant]])
    d <- data.frame(time = time, event = event, score = score)
    d <- d[
      is.finite(d$time) & d$time > 0 & d$event %in% c(0, 1) &
        is.finite(d$score),
    ]
    fit <- coxph(Surv(time, event) ~ score, data = d)
    sm <- summary(fit)
    loo_rows[[length(loo_rows) + 1]] <- data.frame(
      Variant = variant,
      Genes = paste(variants[[variant]], collapse = ";"),
      Cohort = co$name,
      N = nrow(d),
      Events = sum(d$event),
      HR = unname(sm$conf.int["score", "exp(coef)"]),
      CI_low = unname(sm$conf.int["score", "lower .95"]),
      CI_high = unname(sm$conf.int["score", "upper .95"]),
      P = unname(sm$coefficients["score", "Pr(>|z|)"])
    )
  }
}
loo <- do.call(rbind, loo_rows)
write.csv(
  loo,
  file.path(out_dir, "five_gene_leave_one_out_by_cohort.csv"),
  row.names = FALSE
)

meta_fixed <- function(d) {
  beta <- log(d$HR)
  se <- (log(d$CI_high) - log(d$CI_low)) / (2 * 1.96)
  w <- 1 / se^2
  pooled <- sum(w * beta) / sum(w)
  pooled_se <- sqrt(1 / sum(w))
  q <- sum(w * (beta - pooled)^2)
  df <- length(beta) - 1
  i2 <- if (q > 0) max(0, (q - df) / q) * 100 else 0
  data.frame(
    External_meta_HR = exp(pooled),
    CI_low = exp(pooled - 1.96 * pooled_se),
    CI_high = exp(pooled + 1.96 * pooled_se),
    P = 2 * pnorm(-abs(pooled / pooled_se)),
    I2_percent = i2,
    Min_external_HR = min(d$HR),
    Max_external_P = max(d$P)
  )
}

external_names <- c("Moffitt_GSE71729", "Puleo_EMTAB6134", "GSE57495")
loo_meta <- do.call(rbind, lapply(names(variants), function(variant) {
  d <- loo[loo$Variant == variant & loo$Cohort %in% external_names, ]
  cbind(Variant = variant, meta_fixed(d))
}))
loo_meta$P_BH <- p.adjust(loo_meta$P, "BH")
write.csv(
  loo_meta,
  file.path(out_dir, "five_gene_leave_one_out_external_meta.csv"),
  row.names = FALSE
)

gse_probe_sensitivity <- read.csv(
  file.path(out_dir, "GSE57495_probe_aggregation_sensitivity.csv"),
  stringsAsFactors = FALSE
)

report <- c(
  "# 冻结5基因程序：临床多变量与留一基因稳健性",
  paste0("Date: ", Sys.Date()),
  "",
  "Frozen genes: PLOD2, MSN, TPM4, PPFIBP1, SEC31A",
  "",
  "## 临床多变量",
  "- TCGA: age + AJCC stage + grade",
  "- Moffitt: T + N + margin + sex",
  "- Puleo: T + N + grade + margin + sex",
  "- GSE57495: stage",
  "- Clinical_plus_CAF在上述临床模型上增加CAF12等权评分。",
  paste(capture.output(print(multivariable, row.names = FALSE)), collapse = "\n"),
  "",
  "## 留一基因：三外部队列固定效应汇总",
  paste(capture.output(print(loo_meta, row.names = FALSE)), collapse = "\n"),
  "",
  "## GSE57495探针聚合敏感性",
  paste(
    capture.output(print(gse_probe_sensitivity, row.names = FALSE)),
    collapse = "\n"
  ),
  "",
  "## 预设判读",
  "- 若所有Drop-one变体在三个外部队列HR均>1，且外部meta显著，则不是单个基因驱动。",
  "- 不按最优Drop-one结果删基因；主分析仍为冻结Full5。"
)
writeLines(report, file.path(out_dir, "FIVE_GENE_LOCKED_VALIDATION_CN.md"))

print(multivariable)
print(loo_meta)
