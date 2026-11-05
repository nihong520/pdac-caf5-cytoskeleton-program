#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
set.seed(20260820)
suppressPackageStartupMessages(library(WGCNA))
allowWGCNAThreads()
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
pdacr_dir <- file.path(project_root, "data", "raw", "pdacR")
out_dir <- file.path(project_root, "results", "step14_caf_wgcna", "wgcna_standard_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

caf_genes <- c("FAP","ACTA2","COL1A1","COL1A2","COL3A1","PDPN","THY1","DCN","LUM","TAGLN","PDGFRB","SPARC")

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
z_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}
module_score_genes <- function(ex, genes) {
  present <- intersect(toupper(genes), rownames(ex))
  if (!length(present)) return(rep(NA_real_, ncol(ex)))
  as.numeric(colMeans(z_rows(ex[present, , drop = FALSE])))
}
parse_ajcc_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  bad <- x %in% c("", "NA", "NAN", "NULL") |
    grepl("NOT AVAILABLE|DISCREPANCY|UNKNOWN|NOT REPORTED|^\\[", x)
  x[bad] <- NA_character_
  out <- rep(NA_real_, length(x))
  out[!is.na(x) & grepl("IV", x)] <- 4
  out[is.na(out) & !is.na(x) & grepl("III", x)] <- 3
  out[is.na(out) & !is.na(x) & grepl("II", x)] <- 2
  out[is.na(out) & !is.na(x) & grepl("\\bI\\b", x) & !grepl("II|III|IV", x)] <- 1
  out
}
png_hi <- function(file, width = 10, height = 7, res = 300) {
  tryCatch(
    png(file, width = width, height = height, units = "in", res = res, type = "cairo"),
    error = function(e) png(file, width = width, height = height, units = "in", res = res)
  )
}

message("Loading TCGA...")
obj <- readRDS(file.path(pdacr_dir, "TCGA_PAAD.rds"))
keep <- which(
  as.character(obj$sampInfo$Filter_by) == "whitelist" &
    as.character(obj$sampInfo$Histological.type.by.DCC) ==
      "pancreas-adenocarcinoma ductal type"
)
ex <- collapse_symbols(obj, keep)
si <- as.data.frame(obj$sampInfo[keep, , drop = FALSE])
rownames(si) <- colnames(ex)
caf0 <- module_score_genes(ex, caf_genes)
names(caf0) <- colnames(ex)
vars <- apply(ex, 1, var, na.rm = TRUE)
vars[!is.finite(vars)] <- 0
keep_genes <- names(sort(vars, decreasing = TRUE))[seq_len(min(5000L, sum(vars > 0)))]
keep_genes <- union(keep_genes, intersect(caf_genes, rownames(ex)))
datExpr0 <- as.data.frame(t(ex[keep_genes, , drop = FALSE]))
gsg <- goodSamplesGenes(datExpr0, verbose = 0)
if (!gsg$allOK) datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]

sampleTree0 <- hclust(dist(datExpr0), method = "average")
cut_h <- max(sampleTree0$height) * 0.6
clust <- cutreeStatic(sampleTree0, cutHeight = cut_h, minSize = 10)
avg_d <- rowMeans(as.matrix(dist(datExpr0)))
thr_d <- mean(avg_d) + 4 * sd(avg_d)
keep_samp <- (clust == 1) & (avg_d <= thr_d)
removed <- rownames(datExpr0)[!keep_samp]

stage0 <- parse_ajcc_stage(si[rownames(datExpr0), "AJCC.pathologic.tumor.stage"])
traits0 <- data.frame(
  CAF_score = as.numeric(scale(caf0[rownames(datExpr0)])),
  OS_time = suppressWarnings(as.numeric(si[rownames(datExpr0), "survivalA"])),
  OS_event = suppressWarnings(as.numeric(as.character(si[rownames(datExpr0), "censorA.0yes.1no"]))),
  Age = suppressWarnings(as.numeric(as.character(si[rownames(datExpr0), "Age.at.initial.pathologic.diagnosis"]))),
  Stage = stage0,
  Purity = suppressWarnings(as.numeric(as.character(si[rownames(datExpr0), "ABSOLUTE.Purity"]))),
  row.names = rownames(datExpr0)
)
write.csv(
  data.frame(
    Sample = rownames(datExpr0),
    mean_dist = avg_d,
    pass_cutree = (clust == 1),
    pass_dist4sd = (avg_d <= thr_d),
    keep = keep_samp,
    Tumor.Sample.ID = as.character(si[rownames(datExpr0), "Tumor.Sample.ID"])
  ),
  file.path(out_dir, "QC_sample_outlier_flags.csv"),
  row.names = FALSE
)
write.csv(data.frame(Sample = removed), file.path(out_dir, "QC_removed_outliers.csv"), row.names = FALSE)

png_hi(file.path(out_dir, "Fig1A_sample_dendrogram_BEFORE_QC.png"), 11, 7)
plotDendroAndColors(
  sampleTree0, numbers2colors(traits0, signed = FALSE, naColor = "grey80"),
  groupLabels = colnames(traits0),
  main = paste0("Fig1A. BEFORE QC (n=", nrow(datExpr0), "; remove=", sum(!keep_samp), ")"),
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, cex.colorLabels = 0.8
)
abline(h = cut_h, col = "red", lty = 2, lwd = 1.5)
dev.off()

datExpr <- datExpr0[keep_samp, , drop = FALSE]
si <- si[rownames(datExpr), , drop = FALSE]
caf <- caf0[rownames(datExpr)]
message("After QC n=", nrow(datExpr), " removed=", length(removed), ": ", paste(removed, collapse = ", "))

traits <- data.frame(
  CAF_score = as.numeric(scale(caf)),
  OS_time = suppressWarnings(as.numeric(si$survivalA)),
  OS_event = suppressWarnings(as.numeric(as.character(si[["censorA.0yes.1no"]]))),
  Age = suppressWarnings(as.numeric(as.character(si[["Age.at.initial.pathologic.diagnosis"]]))),
  Stage = parse_ajcc_stage(si[["AJCC.pathologic.tumor.stage"]]),
  Purity = suppressWarnings(as.numeric(as.character(si[["ABSOLUTE.Purity"]]))),
  row.names = rownames(datExpr)
)
write.csv(
  as.data.frame(table(Stage = traits$Stage, useNA = "ifany")),
  file.path(out_dir, "QC_stage_encoding_counts.csv"),
  row.names = FALSE
)

sampleTree <- hclust(dist(datExpr), method = "average")
png_hi(file.path(out_dir, "Fig1_sample_dendrogram.png"), 11, 7)
plotDendroAndColors(
  sampleTree, numbers2colors(traits, signed = FALSE, naColor = "grey80"),
  groupLabels = colnames(traits),
  main = paste0("Fig1. AFTER QC (n=", nrow(datExpr), ")"),
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, cex.colorLabels = 0.85
)
dev.off()

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0, networkType = "signed")
fit <- sft$fitIndices
auto_ok <- which(fit$SFT.R.sq >= 0.80)
auto_power <- if (length(auto_ok)) fit$Power[auto_ok[1]] else fit$Power[which.max(fit$SFT.R.sq)]
power <- 14L
message(
  "Auto power=", auto_power, " ; using fixed power=", power,
  " R2=", round(fit$SFT.R.sq[fit$Power == power], 3)
)

png_hi(file.path(out_dir, "Fig2_soft_threshold.png"), 10, 5)
par(mfrow = c(1, 2), mar = c(5, 5, 3.5, 1), cex.lab = 1.1)
plot(fit$Power, fit$SFT.R.sq, xlab = "Soft Threshold (power)", ylab = "Scale-free fit (signed R2)",
     type = "n", main = "A. Scale independence")
text(fit$Power, fit$SFT.R.sq, labels = fit$Power, cex = 0.95, col = "firebrick")
abline(h = 0.80, col = "dodgerblue", lty = 2, lwd = 1.5)
points(power, fit$SFT.R.sq[fit$Power == power], pch = 21, bg = "forestgreen", cex = 1.8)
plot(fit$Power, fit$mean.k., xlab = "Soft Threshold (power)", ylab = "Mean connectivity",
     type = "n", main = "B. Mean connectivity")
text(fit$Power, fit$mean.k., labels = fit$Power, cex = 0.95, col = "firebrick")
points(power, fit$mean.k.[fit$Power == power], pch = 21, bg = "forestgreen", cex = 1.8)
dev.off()

message("Building network power=", power)
net <- blockwiseModules(
  datExpr, power = power, networkType = "signed", TOMType = "signed",
  minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 2
)
moduleColors <- labels2colors(net$colors)
MEs <- orderMEs(net$MEs)
geneTree <- net$dendrograms[[1]]
nSamples <- nrow(datExpr)

moduleTraitCor <- cor(MEs, traits, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
me_ids <- as.integer(sub("^ME", "", colnames(MEs)))
me_colors <- ifelse(me_ids == 0, "grey", labels2colors(me_ids))
rownames(moduleTraitCor) <- paste0(me_colors, " (", colnames(MEs), ")")
rownames(moduleTraitPvalue) <- rownames(moduleTraitCor)

caf_vec <- moduleTraitCor[, "CAF_score"]
non_grey <- !grepl("^grey", rownames(moduleTraitCor))
best_i <- which(non_grey)[which.max(abs(caf_vec[non_grey]))]
best_lab <- rownames(moduleTraitCor)[best_i]
best_color <- me_colors[best_i]
best_num <- me_ids[best_i]
best_r <- unname(caf_vec[best_i])
best_p <- unname(moduleTraitPvalue[best_i, "CAF_score"])
best_n <- sum(moduleColors == best_color)
message(sprintf("Key module=%s n=%d CAF r=%.3f P=%.2e", best_lab, best_n, best_r, best_p))

png_hi(file.path(out_dir, "Fig3_gene_dendrogram_modules.png"), 11, 7)
plotDendroAndColors(
  geneTree, moduleColors[net$blockGenes[[1]]], "Module colors",
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
  main = paste0("Fig3. Gene dendrogram & modules (key=", best_color, " n=", best_n, ")"),
  cex.colorLabels = 0.9
)
dev.off()

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
png_hi(file.path(out_dir, "Fig4_module_trait_heatmap.png"), 9.5, 10)
par(mar = c(7, 11, 3, 2))
labeledHeatmap(
  Matrix = moduleTraitCor, xLabels = colnames(traits),
  yLabels = rownames(moduleTraitCor), ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix,
  setStdMargins = FALSE, cex.text = 0.65, cex.lab = 0.95, zlim = c(-1, 1),
  main = paste0("Fig4. Module-trait (TCGA n=", nSamples, " after QC; key=", best_color, ")")
)
dev.off()
write.csv(cbind(Module = rownames(moduleTraitCor), as.data.frame(moduleTraitCor)),
          file.path(out_dir, "module_trait_cor.csv"), row.names = FALSE)
write.csv(cbind(Module = rownames(moduleTraitPvalue), as.data.frame(moduleTraitPvalue)),
          file.path(out_dir, "module_trait_pvalue.csv"), row.names = FALSE)

modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
names(geneModuleMembership) <- paste0("MM", modNames)
geneTraitSignificance <- as.data.frame(cor(datExpr, traits$CAF_score, use = "p"))
names(geneTraitSignificance) <- "GS.CAF"
column <- match(as.character(best_num), modNames)
moduleGenes <- moduleColors == best_color
mm <- abs(geneModuleMembership[moduleGenes, column])
gs <- abs(geneTraitSignificance[moduleGenes, 1])
gnames <- colnames(datExpr)[moduleGenes]

png_hi(file.path(out_dir, "Fig5_MM_vs_GS_green_CAF.png"), 7.5, 7)
verboseScatterplot(
  mm, gs,
  xlab = paste0("Module Membership |MM| in ", best_color),
  ylab = "Gene significance |GS| for CAF score",
  main = paste0("Fig5. MM vs GS (", best_color, ", n=", sum(moduleGenes), ")"),
  cex.main = 1.15, cex.lab = 1.15, cex.axis = 1.05, col = "darkgreen"
)
abline(v = 0.8, lty = 2, col = "grey40", lwd = 1.2)
abline(h = 0.3, lty = 2, col = "grey40", lwd = 1.2)
hub_idx <- which(mm > 0.8 & gs > 0.3)
if (length(hub_idx)) {
  ord <- order(-(mm[hub_idx] * gs[hub_idx]))
  top <- head(hub_idx[ord], 12)
  text(mm[top], gs[top], labels = gnames[top], cex = 0.6, pos = 3)
}
dev.off()

hub_tab <- data.frame(
  Gene = gnames,
  MM = as.numeric(geneModuleMembership[moduleGenes, column]),
  abs_MM = mm,
  GS_CAF = as.numeric(geneTraitSignificance[moduleGenes, 1]),
  abs_GS = gs,
  ModuleColor = best_color,
  stringsAsFactors = FALSE
)
hub_tab <- hub_tab[order(-hub_tab$abs_MM * hub_tab$abs_GS), ]
hub_tab$Hub_MM08_GS03 <- hub_tab$abs_MM > 0.8 & hub_tab$abs_GS > 0.3
write.csv(hub_tab, file.path(out_dir, "green_MM_GS_all_genes.csv"), row.names = FALSE)
write.csv(hub_tab[hub_tab$Hub_MM08_GS03, ], file.path(out_dir, "green_hub_MM08_GS03.csv"), row.names = FALSE)
write.csv(
  data.frame(Gene = gnames, ModuleColor = best_color),
  file.path(out_dir, "selected_module_genes_AFTER_QC.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    removed = removed, n_before = nrow(datExpr0), n_after = nrow(datExpr),
    power = power, auto_power = auto_power, key_module = best_lab, key_color = best_color,
    key_n = best_n, caf_r = best_r, caf_p = best_p, hub_n = sum(hub_tab$Hub_MM08_GS03)
  ),
  file.path(out_dir, "figure_run_summary_cleaned.rds")
)

lines <- c(
  "# WGCNA figures AFTER QC fixes",
  paste0("Date: ", Sys.Date()),
  "",
  "## Fixes",
  "1. Outlier removal before network (cutreeStatic fail OR mean-dist > mean+4SD).",
  paste0("   Before=", nrow(datExpr0), " After=", nrow(datExpr), " Removed=", length(removed)),
  paste0("   IDs: ", if (length(removed)) paste(removed, collapse = ", ") else "none"),
  "2. AJCC stage encoded I/II/III/IV -> 1/2/3/4.",
  "3. Soft power FIXED at 14 (same as step14 discovery).",
  "4. Key module = max |cor| with CAF among non-grey.",
  "5. 300 dpi redraw; NA traits = grey80.",
  "",
  "## Key results",
  paste0("- Soft power used: ", power, " (auto would be ", auto_power, ")"),
  paste0("- Key module: ", best_lab, " n=", best_n),
  paste0("- CAF r=", signif(best_r, 3), " P=", signif(best_p, 3)),
  paste0("- Hubs |MM|>0.8 & |GS|>0.3: ", sum(hub_tab$Hub_MM08_GS03)),
  "",
  "Use Fig1 (AFTER QC), not Fig1A, in the main figure set."
)
writeLines(lines, file.path(out_dir, "FIGURE_README_CN.md"))
message("Done.")
