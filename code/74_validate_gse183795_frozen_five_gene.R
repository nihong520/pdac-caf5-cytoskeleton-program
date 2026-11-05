#!/usr/bin/env Rscript
# Confirmatory external validation of the FROZEN equal-weight five-gene score
# in GSE183795. Uses submitter gene-level normalized matrix (includes SEC31A).
# Do NOT retune genes/weights based on outcome.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
})

root <- normalizePath(getwd())
out_dir <- file.path(root, "results", "step21_gse183795_frozen_validation")
raw_dir <- file.path(root, "data", "raw", "GSE183795")
pkg_tables <- file.path(
  root, "manuscript", "tables"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pkg_tables, recursive = TRUE, showWarnings = FALSE)

genes <- c("PLOD2", "MSN", "TPM4", "PPFIBP1", "SEC31A")
caf_genes <- c(
  "FAP", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "PDPN",
  "THY1", "DCN", "LUM", "TAGLN", "PDGFRB", "SPARC"
)

norm_gz <- file.path(raw_dir, "GSE183795_normalized_matrix.txt.gz")
soft_gz <- file.path(raw_dir, "GSE183795_family.soft.gz")
stopifnot(file.exists(norm_gz), file.exists(soft_gz))

message("Reading submitter normalized gene matrix...")
norm <- read.delim(gzfile(norm_gz), check.names = FALSE, stringsAsFactors = FALSE)
# Columns: ID, Transcript ID, gene_assignment, RefSeq, then CEL samples
stopifnot("ID" %in% names(norm))
gene_ids <- toupper(trimws(as.character(norm$ID)))
sample_cols <- setdiff(
  names(norm),
  c("ID", "Transcript ID", "gene_assignment", "RefSeq")
)
ex <- as.matrix(norm[, sample_cols, drop = FALSE])
storage.mode(ex) <- "numeric"
rownames(ex) <- gene_ids
# Collapse rare duplicated gene symbols by mean (pre-survival).
if (anyDuplicated(rownames(ex))) {
  message("Collapsing duplicated gene symbols by mean...")
  sums <- rowsum(ex, group = rownames(ex), reorder = FALSE, na.rm = TRUE)
  counts <- as.numeric(table(factor(rownames(ex), levels = rownames(sums))))
  ex <- sums / counts
}
colnames(ex) <- sub("\\.CEL$", "", colnames(ex), ignore.case = TRUE)
message("Gene matrix: ", nrow(ex), " genes x ", ncol(ex), " samples")

coverage <- data.frame(
  Gene = c(genes, caf_genes),
  Set = c(rep("Frozen5", length(genes)), rep("CAF12", length(caf_genes))),
  Present = c(genes, caf_genes) %in% rownames(ex)
)
write.csv(coverage, file.path(out_dir, "GSE183795_gene_coverage.csv"), row.names = FALSE)
write.csv(coverage, file.path(pkg_tables, "GSE183795_gene_coverage.csv"), row.names = FALSE)
print(coverage[coverage$Set == "Frozen5", ])
stopifnot(all(genes %in% rownames(ex)))

message("Parsing phenotype from family soft...")
one_chr <- function(x) {
  if (is.null(x) || length(x) < 1 || is.na(x[1]) || !nzchar(as.character(x[1]))) {
    return(NA_character_)
  }
  as.character(x[1])
}
con <- gzfile(soft_gz, open = "rt")
lines <- readLines(con, warn = FALSE)
close(con)
sample_starts <- grep("^\\^SAMPLE = ", lines)
pheno_list <- vector("list", length(sample_starts))
for (i in seq_along(sample_starts)) {
  start <- sample_starts[i]
  end <- if (i < length(sample_starts)) sample_starts[i + 1] - 1 else length(lines)
  block <- lines[start:end]
  gsm <- one_chr(sub("^\\^SAMPLE = ", "", block[1]))
  get1 <- function(prefix) {
    hit <- grep(paste0("^", prefix, " = "), block, value = TRUE)
    if (!length(hit)) return(NA_character_)
    one_chr(sub(paste0("^", prefix, " = "), "", hit[1]))
  }
  chars <- sub(
    "^!Sample_characteristics_ch1 = ", "",
    grep("^!Sample_characteristics_ch1 = ", block, value = TRUE)
  )
  kv <- list()
  for (ch in chars) {
    if (grepl(":", ch)) {
      key <- tolower(trimws(sub(":.*$", "", ch)))
      kv[[key]] <- one_chr(trimws(sub("^[^:]+:\\s*", "", ch)))
    }
  }
  pheno_list[[i]] <- data.frame(
    Sample = gsm,
    Title = get1("!Sample_title"),
    Tissue = one_chr(kv[["tissue"]]),
    Stage = one_chr(kv[["stage"]]),
    Grade = {
      g <- one_chr(kv[["grading"]])
      if (is.na(g)) one_chr(kv[["grade"]]) else g
    },
    OS_months = suppressWarnings(as.numeric(one_chr(kv[["survival months"]]))),
    Event = suppressWarnings(as.integer(one_chr(kv[["survival status"]]))),
    stringsAsFactors = FALSE
  )
}
pheno <- do.call(rbind, pheno_list)
pheno <- pheno[!is.na(pheno$Title) & nzchar(pheno$Title), , drop = FALSE]
# Soft titles should match CEL basename without .CEL
pheno$key <- pheno$Title
# Prefer unique keys; if duplicate titles, keep first with OS.
pheno <- pheno[order(!is.finite(pheno$OS_months), pheno$Sample), ]
pheno <- pheno[!duplicated(pheno$key), ]
rownames(pheno) <- pheno$key

common <- intersect(colnames(ex), pheno$key)
message(
  "Matched samples: ", length(common),
  " / expr ", ncol(ex), " / pheno ", nrow(pheno)
)
stopifnot(length(common) > 50)
ex <- ex[, common, drop = FALSE]
ph <- pheno[common, , drop = FALSE]

is_tumor <- grepl("^tumor$", ph$Tissue, ignore.case = TRUE)
is_non_tumor <- grepl("adjacent|non[- ]?tumor|normal", ph$Tissue, ignore.case = TRUE)
keep <- is_tumor & !is_non_tumor
if (sum(keep) < 50) {
  keep <- grepl("(^|_).*T$", ph$Title) | grepl("Tp|_T\\.|_T$", ph$Title)
  keep <- keep & !grepl("_N$|Tn|normal|adjacent", ph$Title, ignore.case = TRUE)
  keep <- keep & !grepl("normal|adjacent", ph$Tissue, ignore.case = TRUE)
}
message("Tumor samples: ", sum(keep))

ex_t <- ex[, keep, drop = FALSE]
ph_t <- ph[keep, , drop = FALSE]

score_set <- function(mat, gset) {
  present <- intersect(gset, rownames(mat))
  stopifnot(length(present) == length(gset))
  z <- t(scale(t(mat[present, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  colMeans(z)
}

score <- as.numeric(scale(score_set(ex_t, genes)))
caf_present <- intersect(caf_genes, rownames(ex_t))
caf_score <- if (length(caf_present) >= 8) {
  as.numeric(scale(score_set(ex_t, caf_genes)))
} else {
  rep(NA_real_, ncol(ex_t))
}

roman_map <- c(
  I = 1, IA = 1, IB = 1, II = 2, IIA = 2, IIB = 2,
  III = 3, IIIA = 3, IIIB = 3, IV = 4
)
stage_txt <- toupper(gsub("[^A-Za-z0-9]", "", ph_t$Stage))
stage_num <- suppressWarnings(as.numeric(stage_txt))
stage_from_roman <- unname(roman_map[stage_txt])
stage_final <- ifelse(is.finite(stage_num), stage_num, stage_from_roman)

d <- data.frame(
  Sample = ph_t$Sample,
  Title = ph_t$Title,
  Tissue = ph_t$Tissue,
  OS_months = ph_t$OS_months,
  Event = ph_t$Event,
  Stage = stage_final,
  Grade = ph_t$Grade,
  Score5 = score,
  CAF = caf_score,
  stringsAsFactors = FALSE
)

diag <- data.frame(
  n_expr_cols = ncol(ex),
  n_tumor = nrow(d),
  n_complete_os = sum(
    is.finite(d$OS_months) & d$OS_months > 0 & d$Event %in% c(0, 1),
    na.rm = TRUE
  ),
  events = sum(
    d$Event == 1 & is.finite(d$OS_months) & d$OS_months > 0,
    na.rm = TRUE
  )
)
write.csv(diag, file.path(out_dir, "GSE183795_sample_diagnostics.csv"), row.names = FALSE)
print(diag)

d_surv <- d[is.finite(d$OS_months) & d$OS_months > 0 & d$Event %in% c(0, 1), ]
# Cohort-internal re-standardization after eligibility filter.
d_surv$Score5 <- as.numeric(scale(d_surv$Score5))
if (all(is.finite(d_surv$CAF))) d_surv$CAF <- as.numeric(scale(d_surv$CAF))

write.csv(d_surv, file.path(out_dir, "GSE183795_five_gene_samples.csv"), row.names = FALSE)
write.csv(
  ex_t[genes, match(d_surv$Title, colnames(ex_t)), drop = FALSE],
  file.path(out_dir, "GSE183795_five_gene_expr_matrix.csv")
)

extract_score <- function(fit, model) {
  sm <- summary(fit)
  data.frame(
    Cohort = "GSE183795",
    Model = model,
    N = fit$n,
    Events = fit$nevent,
    HR = unname(sm$conf.int["Score5", "exp(coef)"]),
    CI_low = unname(sm$conf.int["Score5", "lower .95"]),
    CI_high = unname(sm$conf.int["Score5", "upper .95"]),
    P = unname(sm$coefficients["Score5", "Pr(>|z|)"]),
    stringsAsFactors = FALSE
  )
}

fit <- coxph(Surv(OS_months, Event) ~ Score5, data = d_surv)
results <- extract_score(fit, "Univariable")
if (sum(complete.cases(d_surv[, c("Score5", "Stage")])) >= 40) {
  results <- rbind(
    results,
    extract_score(
      coxph(Surv(OS_months, Event) ~ Score5 + Stage, data = d_surv),
      "Stage_adjusted"
    )
  )
}
if (sum(complete.cases(d_surv[, c("Score5", "CAF")])) >= 40) {
  results <- rbind(
    results,
    extract_score(
      coxph(Surv(OS_months, Event) ~ Score5 + CAF, data = d_surv),
      "CAF_adjusted"
    )
  )
}
if (sum(complete.cases(d_surv[, c("Score5", "Stage", "CAF")])) >= 40) {
  results <- rbind(
    results,
    extract_score(
      coxph(Surv(OS_months, Event) ~ Score5 + Stage + CAF, data = d_surv),
      "Stage_CAF_adjusted"
    )
  )
}

d_surv$Risk <- factor(
  ifelse(d_surv$Score5 >= median(d_surv$Score5), "High", "Low"),
  levels = c("Low", "High")
)
lr <- survdiff(Surv(OS_months, Event) ~ Risk, data = d_surv)
km_p <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
results <- rbind(
  results,
  data.frame(
    Cohort = "GSE183795",
    Model = "KM_median_logrank",
    N = nrow(d_surv),
    Events = sum(d_surv$Event),
    HR = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_,
    P = km_p
  )
)

write.csv(results, file.path(out_dir, "GSE183795_five_gene_cox.csv"), row.names = FALSE)
write.csv(results, file.path(pkg_tables, "GSE183795_five_gene_cox.csv"), row.names = FALSE)

km <- survfit(Surv(OS_months, Event) ~ Risk, data = d_surv)
png(
  file.path(out_dir, "GSE183795_KM_median.png"),
  width = 1600, height = 1200, res = 220
)
plot(
  km,
  col = c("#2C7BB6", "#D7191C"),
  lwd = 2,
  xlab = "Overall survival (months)",
  ylab = "Survival probability",
  main = sprintf(
    "GSE183795 frozen equal-weight five-gene score\n(n=%d; log-rank P=%.3g)",
    nrow(d_surv), km_p
  )
)
legend(
  "topright",
  legend = c("Low (below median)", "High (at/above median)"),
  col = c("#2C7BB6", "#D7191C"),
  lwd = 2,
  bty = "n"
)
dev.off()

message("=== GSE183795 frozen five-gene results ===")
print(results)
message("Wrote: ", out_dir)
