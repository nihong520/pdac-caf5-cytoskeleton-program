#!/usr/bin/env Rscript
# Compare the biology-annotated five-gene block with the corresponding
# unsupervised Walktrap community. GSE57495 is deliberately excluded so this
# sensitivity check does not use the final held-out validation cohort.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(survival))

root <- normalizePath(getwd())
pdacr_dir <- file.path(root, "data", "raw", "pdacR")
out_dir <- file.path(root, "results", "step18_yellow317_mcode")

scores <- list(
  Biology5_PLOD2 = c("PLOD2", "MSN", "TPM4", "PPFIBP1", "SEC31A"),
  Walktrap5_VCL = c("VCL", "MSN", "TPM4", "PPFIBP1", "SEC31A")
)
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

module_score <- function(ex, genes) {
  stopifnot(all(genes %in% rownames(ex)))
  z <- t(scale(t(ex[genes, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z)))
}

prepare_cohort <- function(file, name, primary_only = FALSE, tcga = FALSE) {
  obj <- readRDS(file.path(pdacr_dir, file))
  keep <- seq_len(ncol(obj$ex))
  if (primary_only) {
    keep <- which(tolower(as.character(obj$sampInfo$specimen_type)) == "primary")
  }
  if (tcga) {
    keep <- which(
      as.character(obj$sampInfo$Filter_by) == "whitelist" &
        as.character(obj$sampInfo$Histological.type.by.DCC) ==
          "pancreas-adenocarcinoma ductal type"
    )
    oid <- "TCGA-2J-AABV-01A"
    if (oid %in% colnames(obj$ex)) {
      keep <- setdiff(keep, which(colnames(obj$ex) == oid))
    }
  }
  list(
    ex = collapse_symbols(obj, keep),
    si = as.data.frame(obj$sampInfo[keep, , drop = FALSE]),
    name = name
  )
}

cohorts <- list(
  prepare_cohort("TCGA_PAAD.rds", "TCGA_PAAD", tcga = TRUE),
  prepare_cohort(
    "Moffitt_GEO_array.rds", "Moffitt_GSE71729", primary_only = TRUE
  ),
  prepare_cohort("Puleo_array.rds", "Puleo_EMTAB6134")
)

eval_one <- function(co, score_name, genes) {
  score <- module_score(co$ex, genes)
  caf <- module_score(co$ex, intersect(caf_genes, rownames(co$ex)))
  d <- data.frame(
    time = suppressWarnings(as.numeric(co$si$survivalA)),
    event = suppressWarnings(
      as.numeric(as.character(co$si[["censorA.0yes.1no"]]))
    ),
    score = score,
    caf = caf
  )
  d <- d[
    is.finite(d$time) & d$time > 0 & d$event %in% c(0, 1) &
      is.finite(d$score) & is.finite(d$caf),
  ]
  fit <- summary(coxph(Surv(time, event) ~ score, data = d))
  fit_caf <- summary(coxph(Surv(time, event) ~ score + caf, data = d))
  data.frame(
    Score = score_name,
    Cohort = co$name,
    N = nrow(d),
    Events = sum(d$event),
    HR = unname(fit$conf.int["score", "exp(coef)"]),
    CI_low = unname(fit$conf.int["score", "lower .95"]),
    CI_high = unname(fit$conf.int["score", "upper .95"]),
    P = unname(fit$coefficients["score", "Pr(>|z|)"]),
    CAF_HR = unname(fit_caf$conf.int["score", "exp(coef)"]),
    CAF_CI_low = unname(fit_caf$conf.int["score", "lower .95"]),
    CAF_CI_high = unname(fit_caf$conf.int["score", "upper .95"]),
    CAF_P = unname(fit_caf$coefficients["score", "Pr(>|z|)"])
  )
}

res <- do.call(rbind, lapply(names(scores), function(nm) {
  do.call(rbind, lapply(cohorts, function(co) {
    eval_one(co, nm, scores[[nm]])
  }))
}))

meta_one <- function(d, hr, low, high, model) {
  beta <- log(d[[hr]])
  se <- (log(d[[high]]) - log(d[[low]])) / (2 * 1.96)
  w <- 1 / se^2
  pooled <- sum(w * beta) / sum(w)
  pooled_se <- sqrt(1 / sum(w))
  data.frame(
    Score = unique(d$Score),
    Model = model,
    HR = exp(pooled),
    CI_low = exp(pooled - 1.96 * pooled_se),
    CI_high = exp(pooled + 1.96 * pooled_se),
    P = 2 * pnorm(-abs(pooled / pooled_se))
  )
}

external <- res[res$Cohort != "TCGA_PAAD", ]
meta <- do.call(rbind, lapply(names(scores), function(nm) {
  d <- external[external$Score == nm, ]
  rbind(
    meta_one(d, "HR", "CI_low", "CI_high", "Univariable"),
    meta_one(d, "CAF_HR", "CAF_CI_low", "CAF_CI_high", "CAF_adjusted")
  )
}))

write.csv(
  res,
  file.path(out_dir, "biology5_vs_walktrap5_by_cohort.csv"),
  row.names = FALSE
)
write.csv(
  meta,
  file.path(out_dir, "biology5_vs_walktrap5_external_meta.csv"),
  row.names = FALSE
)
print(res)
print(meta)
