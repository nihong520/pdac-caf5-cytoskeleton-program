#!/usr/bin/env Rscript
# Exploratory one-shot prognosis comparison of four pre-defined MCODE1 blocks.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(survival))

project_root <- normalizePath(getwd())
pdacr_dir <- file.path(project_root, "data", "raw", "pdacR")
out_dir <- file.path(project_root, "results", "step18_yellow317_mcode")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

blocks <- list(
  Collagen_scaffold = c(
    "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A2",
    "COL6A3", "COL8A1", "COL10A1", "COL11A1"
  ),
  Basement_ECM = c(
    "FN1", "HSPG2", "LAMA4", "LAMB1", "LAMC1", "NID1", "NID2"
  ),
  Integrin_adhesion = c("ITGAV", "ITGB1", "ITGB5", "PARVA", "VCL"),
  Remodeling_cytoskeleton = c("PLOD2", "MSN", "TPM4", "PPFIBP1", "SEC31A")
)
caf_genes <- c(
  "FAP", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "PDPN",
  "THY1", "DCN", "LUM", "TAGLN", "PDGFRB", "SPARC"
)

write.csv(
  do.call(rbind, lapply(names(blocks), function(nm) {
    data.frame(Block = nm, Gene = blocks[[nm]])
  })),
  file.path(out_dir, "mcode1_four_blocks_locked.csv"),
  row.names = FALSE
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

module_score <- function(ex, genes) {
  present <- intersect(toupper(genes), rownames(ex))
  z <- t(scale(t(ex[present, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  list(score = colMeans(z), genes = present)
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
  list(ex = ex, si = si, name = "TCGA_PAAD")
}

prepare_moffitt <- function() {
  obj <- readRDS(file.path(pdacr_dir, "Moffitt_GEO_array.rds"))
  keep <- which(tolower(as.character(obj$sampInfo$specimen_type)) == "primary")
  si <- as.data.frame(obj$sampInfo[keep, , drop = FALSE])
  ex <- collapse_symbols(obj, keep)
  rownames(si) <- colnames(ex)
  list(ex = ex, si = si, name = "Moffitt_GSE71729")
}

prepare_puleo <- function() {
  obj <- readRDS(file.path(pdacr_dir, "Puleo_array.rds"))
  si <- as.data.frame(obj$sampInfo)
  ex <- collapse_symbols(obj, seq_len(ncol(obj$ex)))
  rownames(si) <- colnames(ex)
  list(ex = ex, si = si, name = "Puleo_EMTAB6134")
}

eval_one <- function(co, block_name, genes) {
  sc <- module_score(co$ex, genes)
  caf <- module_score(co$ex, caf_genes)
  time <- suppressWarnings(as.numeric(co$si$survivalA))
  event <- suppressWarnings(as.numeric(as.character(co$si[["censorA.0yes.1no"]])))
  purity <- if ("ABSOLUTE.Purity" %in% names(co$si)) {
    suppressWarnings(as.numeric(as.character(co$si$ABSOLUTE.Purity)))
  } else {
    rep(NA_real_, nrow(co$si))
  }
  d <- data.frame(
    time = time,
    event = event,
    score = as.numeric(scale(sc$score)),
    caf = as.numeric(scale(caf$score)),
    purity = purity
  )
  d <- d[
    is.finite(d$time) & d$time > 0 & d$event %in% c(0, 1) &
      is.finite(d$score) & is.finite(d$caf),
  ]

  fit <- coxph(Surv(time, event) ~ score, data = d)
  sm <- summary(fit)
  fit_caf <- coxph(Surv(time, event) ~ score + caf, data = d)
  sm_caf <- summary(fit_caf)

  purity_hr <- purity_p <- NA_real_
  if (sum(is.finite(d$purity)) >= 40) {
    dp <- d[is.finite(d$purity), ]
    fit_purity <- coxph(Surv(time, event) ~ score + purity, data = dp)
    sm_purity <- summary(fit_purity)
    purity_hr <- unname(sm_purity$conf.int["score", "exp(coef)"])
    purity_p <- unname(sm_purity$coefficients["score", "Pr(>|z|)"])
  }

  data.frame(
    Block = block_name,
    Cohort = co$name,
    N = nrow(d),
    Events = sum(d$event),
    Genes_used = length(sc$genes),
    HR = unname(sm$conf.int["score", "exp(coef)"]),
    CI_low = unname(sm$conf.int["score", "lower .95"]),
    CI_high = unname(sm$conf.int["score", "upper .95"]),
    P = unname(sm$coefficients["score", "Pr(>|z|)"]),
    CAF_adj_HR = unname(sm_caf$conf.int["score", "exp(coef)"]),
    CAF_adj_CI_low = unname(sm_caf$conf.int["score", "lower .95"]),
    CAF_adj_CI_high = unname(sm_caf$conf.int["score", "upper .95"]),
    CAF_adj_P = unname(sm_caf$coefficients["score", "Pr(>|z|)"]),
    Purity_adj_HR = purity_hr,
    Purity_adj_P = purity_p,
    Rho_vs_CAF = cor(d$score, d$caf, use = "complete.obs")
  )
}

message("Loading cohorts...")
cohorts <- list(prepare_tcga(), prepare_moffitt(), prepare_puleo())
res <- do.call(rbind, lapply(names(blocks), function(block_name) {
  do.call(rbind, lapply(cohorts, function(co) {
    eval_one(co, block_name, blocks[[block_name]])
  }))
}))

# Four block tests are adjusted within each cohort/model family.
res$P_BH <- ave(res$P, res$Cohort, FUN = function(x) p.adjust(x, "BH"))
res$CAF_adj_P_BH <- ave(
  res$CAF_adj_P, res$Cohort, FUN = function(x) p.adjust(x, "BH")
)

# Fixed-effect meta-analysis across the two external cohorts.
external <- res[res$Cohort != "TCGA_PAAD", ]
meta_one <- function(d, hr_col, low_col, high_col, model) {
  beta <- log(d[[hr_col]])
  se <- (log(d[[high_col]]) - log(d[[low_col]])) / (2 * 1.96)
  w <- 1 / se^2
  pooled_beta <- sum(w * beta) / sum(w)
  pooled_se <- sqrt(1 / sum(w))
  z <- pooled_beta / pooled_se
  q <- sum(w * (beta - pooled_beta)^2)
  df <- length(beta) - 1
  i2 <- if (q > 0) max(0, (q - df) / q) * 100 else 0
  data.frame(
    Block = unique(d$Block),
    Model = model,
    External_meta_HR = exp(pooled_beta),
    CI_low = exp(pooled_beta - 1.96 * pooled_se),
    CI_high = exp(pooled_beta + 1.96 * pooled_se),
    P = 2 * pnorm(-abs(z)),
    I2_percent = i2
  )
}
meta <- do.call(rbind, lapply(names(blocks), function(block_name) {
  d <- external[external$Block == block_name, ]
  rbind(
    meta_one(d, "HR", "CI_low", "CI_high", "Univariable"),
    meta_one(
      d, "CAF_adj_HR", "CAF_adj_CI_low", "CAF_adj_CI_high",
      "CAF_adjusted"
    )
  )
}))
meta$P_BH <- ave(meta$P, meta$Model, FUN = function(x) p.adjust(x, "BH"))

write.csv(
  res,
  file.path(out_dir, "mcode1_four_blocks_cox.csv"),
  row.names = FALSE
)
write.csv(
  meta,
  file.path(out_dir, "mcode1_four_blocks_external_meta.csv"),
  row.names = FALSE
)

report <- c(
  "# MCODE1四个功能分块的探索性预后分析",
  paste0("Date: ", Sys.Date()),
  "",
  "## 设计",
  "- 四组按生物学注释预先冻结，不按生存结果改组。",
  "- 每组采用基因逐行z标准化后的等权均值。",
  "- 三队列单因素Cox；CAF校正；TCGA另做纯度校正。",
  "- 四组P值在每个队列内使用BH校正。",
  "- Moffitt与Puleo采用固定效应作外部汇总，仅作方向性旁证。",
  "",
  "## 三队列结果",
  paste(capture.output(print(res, row.names = FALSE)), collapse = "\n"),
  "",
  "## 两个外部队列固定效应汇总",
  paste(capture.output(print(meta, row.names = FALSE)), collapse = "\n"),
  "",
  "## 判读规则",
  "- 优先看两个外部队列方向是否一致、外部meta的BH校正及CAF校正。",
  "- 单队列P<0.05但外部不一致或CAF校正消失，不能称独立预后程序。",
  "- 本分析是post-hoc功能分块探索，不用于继续挑基因建模。"
)
writeLines(
  report,
  file.path(out_dir, "MCODE1_FOUR_BLOCKS_PROGNOSIS_CN.md")
)

print(res)
print(meta)
