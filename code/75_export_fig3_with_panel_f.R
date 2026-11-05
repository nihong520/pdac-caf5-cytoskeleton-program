#!/usr/bin/env Rscript
# Rebuild Fig.3: A forest, B leave-one-gene, C leave-one-cohort,
# D equal vs TCGA-weighted, E GSE183795 confirmatory.
# Incremental ΔC is Supplementary Fig. S3.
options(stringsAsFactors = FALSE)
set.seed(20260821)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
if (!is.na(file_arg) && nzchar(file_arg)) {
  script_path <- normalizePath(file_arg, winslash = "/")
  root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
} else {
  root <- normalizePath(getwd(), winslash = "/")
}

pkg <- file.path(root, "manuscript")
main_dir <- file.path(pkg, "figures", "main")
out_dir <- file.path(pkg, "figures", "results_section2_4_ordered")
sub_main <- file.path(pkg, "submission_pack", "figures", "main")
supp_dir <- file.path(pkg, "figures", "supplementary")
sub_supp <- file.path(pkg, "submission_pack", "figures", "supplementary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sub_main, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sub_supp, recursive = TRUE, showWarnings = FALSE)

step19 <- file.path(root, "results", "step19_five_gene_locked_validation")
step20 <- file.path(root, "results", "step20_statistical_strengthening")
step21 <- file.path(root, "results", "step21_gse183795_frozen_validation")
step25 <- file.path(root, "results", "step25_single_gene_prognosis_proxy")

theme_journal <- function(base = 8) {
  theme_classic(base_size = base) %+replace%
    theme(
      plot.title = element_text(face = "bold", size = 8, hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 8, hjust = 0, margin = margin(b = 4)),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "grey94", color = "grey70"),
      plot.tag = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      complete = TRUE
    )
}

normalize_svg_text <- function(svg_path) {
  svg_lines <- readLines(svg_path, warn = FALSE, encoding = "UTF-8")
  is_text <- grepl("<text\\b", svg_lines)
  is_panel <- is_text &
    grepl(">[A-G]</text>", svg_lines) &
    grepl("font-weight: bold|font-weight='bold'|font-weight: 700|font-weight:700", svg_lines)
  svg_lines[is_text & !is_panel] <- gsub(
    "font-size: [0-9.]+px", "font-size: 8.00pt",
    svg_lines[is_text & !is_panel]
  )
  svg_lines[is_panel] <- gsub(
    "font-size: [0-9.]+px", "font-size: 10.00pt",
    svg_lines[is_panel]
  )
  svg_lines[is_text] <- gsub(
    " textLength='[^']*' lengthAdjust='spacingAndGlyphs'",
    "", svg_lines[is_text]
  )
  writeLines(svg_lines, svg_path, useBytes = TRUE)
}

save_vector <- function(plot, stem, width, height, also_main = NULL, also_supp = NULL) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(out_dir, paste0(stem, ".svg"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(
    svg_path, plot, width = width, height = height,
    device = svglite::svglite, bg = "white"
  )
  ggsave(png_path, plot, width = width, height = height, dpi = 300, bg = "white")
  normalize_svg_text(svg_path)
  if (!is.null(also_main)) {
    for (td in unique(c(main_dir, sub_main))) {
      file.copy(pdf_path, file.path(td, paste0(also_main, ".pdf")), overwrite = TRUE)
      file.copy(svg_path, file.path(td, paste0(also_main, ".svg")), overwrite = TRUE)
      file.copy(png_path, file.path(td, paste0(also_main, ".png")), overwrite = TRUE)
    }
  }
  if (!is.null(also_supp)) {
    for (td in unique(c(supp_dir, sub_supp))) {
      file.copy(pdf_path, file.path(td, paste0(also_supp, ".pdf")), overwrite = TRUE)
      file.copy(svg_path, file.path(td, paste0(also_supp, ".svg")), overwrite = TRUE)
      file.copy(png_path, file.path(td, paste0(also_supp, ".png")), overwrite = TRUE)
    }
  }
  message("Wrote ", stem)
}

# ---- Panels A–E (same content as 71_export...) ----
forest <- read.csv(file.path(step19, "figure_forest_three_models_data.csv"), stringsAsFactors = FALSE)
forest$Model <- factor(
  forest$Model,
  levels = c("Univariable", "Clinical adjusted", "Clinical + CAF adjusted")
)
forest$Cohort <- factor(
  forest$Cohort,
  levels = rev(c(
    "TCGA-PAAD", "Moffitt (GSE71729)", "Puleo (E-MTAB-6134)",
    "GSE57495", "External meta"
  ))
)
forest$Role <- factor(forest$Role, levels = c("Discovery", "Validation", "Meta"))
forest$Text <- sprintf(
  "%.2f (%.2f–%.2f)\n%s",
  forest$HR, forest$CI_low, forest$CI_high,
  ifelse(forest$P < 0.001, "P<0.001", sprintf("P=%.3f", forest$P))
)
role_colors <- c(Discovery = "#4D4D4D", Validation = "#0072B2", Meta = "#D55E00")
p_a <- ggplot(forest, aes(x = HR, y = Cohort, color = Role, shape = Role)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0, linewidth = 0.65) +
  geom_point(size = 2.2) +
  geom_text(aes(x = 3.55, label = Text), hjust = 0, color = "black", size = 2.3, lineheight = 0.9) +
  facet_wrap(~Model, nrow = 1) +
  scale_x_log10(limits = c(0.8, 7.2), breaks = c(1, 1.5, 2, 3, 5)) +
  scale_color_manual(values = role_colors) +
  scale_shape_manual(values = c(Discovery = 16, Validation = 16, Meta = 18)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Frozen five-gene score and OS",
    x = "Hazard ratio (log scale)", y = NULL, color = NULL, shape = NULL
  ) +
  theme_journal() +
  theme(plot.margin = margin(5, 70, 5, 5))

loo_gene <- read.csv(file.path(step19, "five_gene_leave_one_out_external_meta.csv"), stringsAsFactors = FALSE)
loo_gene$Label <- c(
  Full5 = "Full five-gene score",
  Drop_PLOD2 = "Without PLOD2", Drop_MSN = "Without MSN",
  Drop_TPM4 = "Without TPM4", Drop_PPFIBP1 = "Without PPFIBP1",
  Drop_SEC31A = "Without SEC31A"
)[loo_gene$Variant]
loo_gene$Label <- factor(loo_gene$Label, levels = rev(loo_gene$Label))
p_b <- ggplot(loo_gene, aes(x = External_meta_HR, y = Label)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0, linewidth = 0.65, color = "#2C7FB8") +
  geom_point(size = 2.2, color = "#2C7FB8") +
  scale_x_log10(limits = c(0.95, 1.8), breaks = c(1, 1.2, 1.4, 1.6)) +
  labs(title = "Leave-one-gene-out", x = "External pooled HR (log)", y = NULL) +
  theme_journal()

loco <- read.csv(file.path(step20, "five_gene_leave_one_cohort_out_meta.csv"), stringsAsFactors = FALSE)
loco$ModelLabel <- ifelse(loco$Model == "Clinical_adjusted", "Clinical", "Clinical + CAF")
loco$Dropped <- c(
  Moffitt_GSE71729 = "Drop Moffitt",
  Puleo_EMTAB6134 = "Drop Puleo",
  GSE57495 = "Drop GSE57495"
)[loco$Dropped_cohort]
loco$Label <- factor(
  paste(loco$ModelLabel, loco$Dropped, sep = " | "),
  levels = rev(c(
    "Clinical | Drop Moffitt", "Clinical | Drop Puleo", "Clinical | Drop GSE57495",
    "Clinical + CAF | Drop Moffitt", "Clinical + CAF | Drop Puleo", "Clinical + CAF | Drop GSE57495"
  ))
)
p_c <- ggplot(loco, aes(x = HR, y = Label, color = ModelLabel)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0, linewidth = 0.65) +
  geom_point(size = 2.1) +
  scale_x_log10(limits = c(0.95, 2.5), breaks = c(1, 1.25, 1.5, 2)) +
  scale_color_manual(values = c(Clinical = "#0072B2", `Clinical + CAF` = "#D55E00")) +
  labs(title = "Leave-one-cohort-out", x = "Pooled HR (log)", y = NULL, color = NULL) +
  theme_journal()

inc <- read.csv(file.path(step20, "five_gene_incremental_value.csv"), stringsAsFactors = FALSE)
inc$CohortLabel <- c(
  TCGA_PAAD = "TCGA", Moffitt_GSE71729 = "Moffitt",
  Puleo_EMTAB6134 = "Puleo", GSE57495 = "GSE57495"
)[inc$Cohort]
inc$Family <- factor(
  inc$Model_family,
  levels = c("Clinical", "Clinical_plus_CAF"),
  labels = c("Clinical", "Clinical + CAF")
)
inc$Label <- factor(
  paste(inc$CohortLabel, inc$Family, sep = " | "),
  levels = rev(as.vector(outer(
    c("TCGA", "Moffitt", "Puleo", "GSE57495"),
    c("Clinical", "Clinical + CAF"), paste, sep = " | "
  )))
)
p_d <- ggplot(inc, aes(x = Delta_C, y = Label, color = Family)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbar(aes(xmin = Delta_C_boot_low, xmax = Delta_C_boot_high), orientation = "y", width = 0, linewidth = 0.6) +
  geom_point(aes(shape = LRT_P < 0.05), size = 2.3) +
  scale_color_manual(values = c(Clinical = "#0072B2", `Clinical + CAF` = "#D55E00")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(
    title = "Incremental C-index",
    x = expression(Delta * " Harrell C-index"), y = NULL, color = NULL
  ) +
  theme_journal()

weight <- read.csv(
  file.path(step19, "equal_vs_TCGAweight_external_meta.csv"),
  stringsAsFactors = FALSE
)
weight$Definition <- factor(
  c(EqualWeight = "Equal weight", TCGAweighted = "TCGA-weighted")[weight$Score],
  levels = c("TCGA-weighted", "Equal weight")
)
weight$Model <- factor(
  weight$Model, levels = c("Univariable", "CAF_adjusted"),
  labels = c("Univariable", "CAF adjusted")
)
p_e <- ggplot(weight, aes(x = HR, y = Definition, color = Model, shape = Model)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60") +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0, linewidth = 0.6,
    position = position_dodge(width = 0.45)
  ) +
  geom_point(position = position_dodge(width = 0.45), size = 2.2) +
  scale_x_log10(limits = c(0.95, 1.9), breaks = c(1, 1.2, 1.5, 1.8)) +
  scale_color_manual(values = c(Univariable = "#0072B2", `CAF adjusted` = "#D55E00")) +
  scale_shape_manual(values = c(Univariable = 16, `CAF adjusted` = 17)) +
  labs(
    title = "Equal vs TCGA-weighted score",
    x = "External pooled HR (log)", y = NULL, color = NULL, shape = NULL
  ) +
  theme_journal()

# ---- Panel F: GSE183795 confirmatory ----
g183 <- read.csv(file.path(step21, "GSE183795_five_gene_cox.csv"), stringsAsFactors = FALSE)
g183 <- g183[g183$Model != "KM_median_logrank", , drop = FALSE]
g183$ModelLabel <- factor(
  c(
    Univariable = "Univariable",
    Stage_adjusted = "Stage adjusted",
    CAF_adjusted = "CAF adjusted",
    Stage_CAF_adjusted = "Stage + CAF"
  )[g183$Model],
  levels = rev(c("Univariable", "Stage adjusted", "CAF adjusted", "Stage + CAF"))
)
g183$Text <- sprintf(
  "%.2f (%.2f–%.2f)\n%s",
  g183$HR, g183$CI_low, g183$CI_high,
  ifelse(g183$P < 0.001, "P<0.001", sprintf("P=%.3f", g183$P))
)
write.csv(
  g183[, c("Cohort", "Model", "N", "Events", "HR", "CI_low", "CI_high", "P")],
  file.path(pkg, "tables", "GSE183795_five_gene_cox.csv"),
  row.names = FALSE
)

p_f <- ggplot(g183, aes(x = HR, y = ModelLabel)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0, linewidth = 0.65, color = "#0072B2") +
  geom_point(size = 2.2, color = "#0072B2") +
  geom_text(aes(x = 2.85, label = Text), hjust = 0, color = "black", size = 2.2, lineheight = 0.9) +
  scale_x_log10(limits = c(0.9, 5.5), breaks = c(1, 1.5, 2, 3)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "GSE183795 confirmatory",
    subtitle = sprintf("n=%d; %d events; frozen equal-weight score", g183$N[1], g183$Events[1]),
    x = "HR per 1 SD (log)", y = NULL
  ) +
  theme_journal() +
  theme(plot.margin = margin(5, 55, 5, 5))

# Main Fig. 3 is A–E: forest, leave-one-gene, leave-one-cohort, equal vs weighted, GSE183795.
# Incremental ΔC is Supplementary Fig. S3.
fig3 <- p_a / (p_b | p_c) / (p_e | p_f) +
  plot_layout(heights = c(1.35, 1, 1.15)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 10, face = "bold"))
  ) &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_vector(
  fig3, "03_Figure3_prognostic_validation_robustness",
  width = 14.5, height = 12.4,
  also_main = "Figure3_prognostic_validation_robustness"
)

message("Fig.3 A–E export complete (ΔC is Supplementary Fig. S3).")

# ---- Supplementary Fig. S3: incremental C-index only ----
cohort_lv <- c("TCGA", "Moffitt", "Puleo", "GSE57495")
fam_cols <- c(Clinical = "#F8766D", `Clinical + CAF` = "#00BFC4")
inc_s3 <- inc
inc_s3$Cohort <- factor(
  c(
    TCGA_PAAD = "TCGA", Moffitt_GSE71729 = "Moffitt",
    Puleo_EMTAB6134 = "Puleo", GSE57495 = "GSE57495"
  )[inc_s3$Cohort],
  levels = rev(cohort_lv)
)
inc_s3$Family <- factor(
  inc_s3$Model_family,
  levels = c("Clinical", "Clinical_plus_CAF"),
  labels = c("Clinical", "Clinical + CAF")
)
fig_s3 <- ggplot(inc_s3, aes(x = Delta_C, y = Cohort, color = Family)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbar(
    aes(xmin = Delta_C_boot_low, xmax = Delta_C_boot_high),
    orientation = "y", width = 0, linewidth = 0.55,
    position = position_dodge(width = 0.45)
  ) +
  geom_point(aes(shape = LRT_P < 0.05), position = position_dodge(width = 0.45), size = 2.2) +
  scale_color_manual(values = fam_cols) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_x_continuous(
    limits = c(-0.05, 0.26),
    breaks = c(-0.05, 0, 0.05, 0.10, 0.15, 0.20, 0.25)
  ) +
  labs(
    title = "Incremental C-index",
    x = expression(Delta * " Harrell C-index"),
    y = NULL, color = NULL,
    caption = paste(
      "Fig. S3 | Change in Harrell's C-index after adding the five-gene score",
      "to clinical or clinical+CAF models. Lines are 1,000-bootstrap percentile 95% CIs;",
      "filled points indicate a likelihood-ratio test P<0.05."
    )
  ) +
  theme_journal() +
  theme(
    legend.position = "bottom",
    plot.caption = element_text(size = 8, hjust = 0, lineheight = 1.2, margin = margin(t = 8))
  )
save_vector(
  fig_s3, "04_FigureS3_model_diagnostics_platform_sensitivity",
  width = 6.2, height = 4.2,
  also_supp = "FigureS3_model_diagnostics_platform_sensitivity"
)
message("Supplementary Fig. S3 export complete.")
