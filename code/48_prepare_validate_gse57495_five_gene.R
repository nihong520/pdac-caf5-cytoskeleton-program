#!/usr/bin/env Rscript
# Independent validation of the frozen PLOD2/MSN/TPM4/PPFIBP1/SEC31A score.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(survival)
})

project_root <- normalizePath(getwd())
gse57495_candidates <- c(
  Sys.getenv("GSE57495_SERIES_MATRIX", unset = ""),
  file.path(project_root, "data", "raw", "GSE57495", "GSE57495_series_matrix.txt")
)
gse57495_candidates <- gse57495_candidates[nzchar(gse57495_candidates)]
matrix_file <- gse57495_candidates[file.exists(gse57495_candidates)][1]
if (is.na(matrix_file) || !nzchar(matrix_file)) {
  stop("Place the GEO GSE57495 series matrix at data/raw/GSE57495/GSE57495_series_matrix.txt")
}
out_dir <- file.path(project_root, "results", "step19_five_gene_locked_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
genes <- c("PLOD2", "MSN", "TPM4", "PPFIBP1", "SEC31A")
caf_genes <- c(
  "FAP", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "PDPN",
  "THY1", "DCN", "LUM", "TAGLN", "PDGFRB", "SPARC"
)

message("Reading local GSE57495 series matrix...")
eset <- getGEO(filename = matrix_file, getGPL = FALSE)
if (is.list(eset) && !methods::is(eset, "ExpressionSet")) eset <- eset[[1]]
probe_expr <- exprs(eset)
pd <- pData(eset)
message("Expression: ", nrow(probe_expr), " probes x ", ncol(probe_expr), " samples")

gpl_data_file <- file.path(out_dir, "GPL15048_data.soft")
gpl_cache <- file.path(out_dir, "GPL15048_annotation.rds")
if (file.exists(gpl_data_file)) {
  message("Reading cached platform-only GPL15048 table")
  gpl_lines <- readLines(gpl_data_file, warn = FALSE)
  table_start <- grep("^!platform_table_begin", gpl_lines)
  table_end <- grep("^!platform_table_end", gpl_lines)
  if (!length(table_start) || !length(table_end)) {
    stop("Cached GPL15048 table is incomplete")
  }
  gpl_tab <- read.delim(
    text = gpl_lines[(table_start + 1):(table_end - 1)],
    check.names = FALSE,
    quote = "",
    stringsAsFactors = FALSE
  )
  saveRDS(gpl_tab, gpl_cache)
} else if (!file.exists(gpl_cache)) {
  message("Downloading GPL15048 annotation from GEO...")
  gpl <- getGEO("GPL15048", destdir = out_dir, AnnotGPL = TRUE)
  gpl_tab <- Table(gpl)
  saveRDS(gpl_tab, gpl_cache)
} else {
  message("Using cached GPL15048 annotation")
  gpl_tab <- readRDS(gpl_cache)
}

id_col <- grep("^ID$|ID_REF|probe", names(gpl_tab), ignore.case = TRUE, value = TRUE)[1]
sym_col <- grep(
  "^Gene Symbol$|gene.?symbol|^symbol$",
  names(gpl_tab), ignore.case = TRUE, value = TRUE
)[1]
if (is.na(id_col) || is.na(sym_col)) {
  stop(
    "GPL annotation lacks ID/SYMBOL columns. Available: ",
    paste(names(gpl_tab), collapse = ", ")
  )
}

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
if (!any(keep)) stop("No probes mapped to gene symbols")

# Select the highest-variance probe per gene, a platform-stable pre-specified rule.
mapped_expr <- probe_expr[keep, , drop = FALSE]
mapped_sym <- mapped[keep]
# Sensitivity matrix: mean all probes mapping to the same gene.
mean_sums <- rowsum(mapped_expr, group = mapped_sym, reorder = FALSE, na.rm = TRUE)
mean_counts <- as.numeric(table(factor(mapped_sym, levels = rownames(mean_sums))))
mean_expr <- mean_sums / mean_counts

# Main matrix: highest-variance probe per gene, selected without survival data.
probe_var <- apply(mapped_expr, 1, var, na.rm = TRUE)
ord <- order(probe_var, decreasing = TRUE, na.last = TRUE)
mapped_expr <- mapped_expr[ord, , drop = FALSE]
mapped_sym <- mapped_sym[ord]
mapped_expr <- mapped_expr[!duplicated(mapped_sym), , drop = FALSE]
rownames(mapped_expr) <- mapped_sym[!duplicated(mapped_sym)]

message("Gene-level matrix: ", nrow(mapped_expr), " genes")
coverage <- data.frame(
  Gene = c(genes, caf_genes),
  Set = c(rep("Frozen5", length(genes)), rep("CAF12", length(caf_genes))),
  Present = c(genes, caf_genes) %in% rownames(mapped_expr)
)
write.csv(coverage, file.path(out_dir, "GSE57495_gene_coverage.csv"), row.names = FALSE)
message("Frozen genes present: ", paste(genes[genes %in% rownames(mapped_expr)], collapse = ", "))

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
stage_raw <- extract_characteristic(pd, "^stage\\s*:")
os_months <- suppressWarnings(as.numeric(sub(".*:\\s*", "", os_raw)))
event <- ifelse(grepl("dead", event_raw, ignore.case = TRUE), 1,
                ifelse(grepl("alive", event_raw, ignore.case = TRUE), 0, NA))
stage_text <- toupper(sub(".*:\\s*", "", stage_raw))
stage_num <- suppressWarnings(as.numeric(sub("^([1-4]).*$", "\\1", stage_text)))

score_set <- function(mat, gset) {
  present <- intersect(gset, rownames(mat))
  z <- t(scale(t(mat[present, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  colMeans(z)
}
score <- as.numeric(scale(score_set(mapped_expr, genes)))
caf_present <- intersect(caf_genes, rownames(mapped_expr))
caf_score <- if (length(caf_present) >= 5) {
  as.numeric(scale(score_set(mapped_expr, caf_genes)))
} else {
  rep(NA_real_, ncol(mapped_expr))
}

d <- data.frame(
  Sample = colnames(mapped_expr),
  OS_months = os_months,
  Event = event,
  Stage = stage_num,
  Score5 = score,
  CAF = caf_score
)
d <- d[is.finite(d$OS_months) & d$OS_months > 0 & d$Event %in% c(0, 1), ]
write.csv(d, file.path(out_dir, "GSE57495_five_gene_samples.csv"), row.names = FALSE)

fit <- coxph(Surv(OS_months, Event) ~ Score5, data = d)
sm <- summary(fit)
extract_score <- function(sm, model) {
  data.frame(
    Cohort = "GSE57495",
    Model = model,
    N = sm$n,
    Events = sm$nevent,
    HR = unname(sm$conf.int["Score5", "exp(coef)"]),
    CI_low = unname(sm$conf.int["Score5", "lower .95"]),
    CI_high = unname(sm$conf.int["Score5", "upper .95"]),
    P = unname(sm$coefficients["Score5", "Pr(>|z|)"])
  )
}
results <- extract_score(sm, "Univariable")

if (sum(complete.cases(d[, c("Score5", "Stage")])) >= 40) {
  fit_stage <- coxph(Surv(OS_months, Event) ~ Score5 + Stage, data = d)
  results <- rbind(results, extract_score(summary(fit_stage), "Stage_adjusted"))
}
if (sum(complete.cases(d[, c("Score5", "CAF")])) >= 40) {
  fit_caf <- coxph(Surv(OS_months, Event) ~ Score5 + CAF, data = d)
  results <- rbind(results, extract_score(summary(fit_caf), "CAF_adjusted"))
}

write.csv(results, file.path(out_dir, "GSE57495_five_gene_cox.csv"), row.names = FALSE)

# Probe-aggregation sensitivity: repeat with the mean of all mapped probes.
mean_score <- as.numeric(scale(score_set(mean_expr, genes)))
mean_caf_score <- if (length(intersect(caf_genes, rownames(mean_expr))) >= 5) {
  as.numeric(scale(score_set(mean_expr, caf_genes)))
} else {
  rep(NA_real_, ncol(mean_expr))
}
dm <- data.frame(
  OS_months = os_months,
  Event = event,
  Stage = stage_num,
  Score5 = mean_score,
  CAF = mean_caf_score
)
dm <- dm[
  is.finite(dm$OS_months) & dm$OS_months > 0 & dm$Event %in% c(0, 1),
]
mean_results <- extract_score(
  summary(coxph(Surv(OS_months, Event) ~ Score5, data = dm)),
  "Mean_probes_univariable"
)
mean_results <- rbind(
  mean_results,
  extract_score(
    summary(coxph(Surv(OS_months, Event) ~ Score5 + Stage, data = dm)),
    "Mean_probes_stage_adjusted"
  )
)
if (sum(complete.cases(dm[, c("Score5", "CAF")])) >= 40) {
  mean_results <- rbind(
    mean_results,
    extract_score(
      summary(coxph(Surv(OS_months, Event) ~ Score5 + CAF, data = dm)),
      "Mean_probes_CAF_adjusted"
    )
  )
}
probe_sensitivity <- rbind(
  transform(results, Aggregation = "Highest_variance_probe"),
  transform(mean_results, Aggregation = "Mean_all_probes")
)
write.csv(
  probe_sensitivity,
  file.path(out_dir, "GSE57495_probe_aggregation_sensitivity.csv"),
  row.names = FALSE
)
saveRDS(mapped_expr, file.path(out_dir, "GSE57495_expr_genes.rds"))
print(results)
print(mean_results)
