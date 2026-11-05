#!/usr/bin/env Rscript
# Publication figures for the locked 5-gene score.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

project_root <- normalizePath(getwd())
step18 <- file.path(project_root, "results", "step18_yellow317_mcode")
out_dir <- file.path(
  project_root, "results", "step19_five_gene_locked_validation"
)
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

genes5 <- c("PLOD2", "MSN", "TPM4", "PPFIBP1", "SEC31A")
write.csv(
  data.frame(
    Gene = genes5,
    Weight = rep(0.2, length(genes5)),
    Score = "remodeling_cytoskeleton_5gene_v1.0",
    Locked_date = "2026-08-21"
  ),
  file.path(out_dir, "FIVE_GENE_SCORE_SPECIFICATION_LOCKED.csv"),
  row.names = FALSE
)

fixed_meta <- function(d) {
  beta <- log(d$HR)
  se <- (log(d$CI_high) - log(d$CI_low)) / (2 * 1.96)
  w <- 1 / se^2
  pooled <- sum(w * beta) / sum(w)
  pooled_se <- sqrt(1 / sum(w))
  q <- sum(w * (beta - pooled)^2)
  df <- length(beta) - 1
  i2 <- if (q > 0) max(0, (q - df) / q) * 100 else 0
  data.frame(
    Cohort = "External meta",
    N = sum(d$N),
    Events = sum(d$Events),
    HR = exp(pooled),
    CI_low = exp(pooled - 1.96 * pooled_se),
    CI_high = exp(pooled + 1.96 * pooled_se),
    P = 2 * pnorm(-abs(pooled / pooled_se)),
    I2_percent = i2
  )
}

format_p <- function(p) {
  ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
         sprintf("%.3f", p))
}

# Assemble univariable, clinical-adjusted, and clinical+CAF models.
u0 <- read.csv(
  file.path(step18, "mcode1_four_blocks_cox.csv"),
  stringsAsFactors = FALSE
)
u0 <- u0[u0$Block == "Remodeling_cytoskeleton", ]
univ <- u0[, c("Cohort", "N", "Events", "HR", "CI_low", "CI_high", "P")]
g57 <- read.csv(
  file.path(out_dir, "GSE57495_five_gene_cox.csv"),
  stringsAsFactors = FALSE
)
g57 <- g57[g57$Model == "Univariable", ]
univ <- rbind(
  univ,
  g57[, c("Cohort", "N", "Events", "HR", "CI_low", "CI_high", "P")]
)
univ$Model <- "Univariable"

multi <- read.csv(
  file.path(out_dir, "five_gene_clinical_multivariable.csv"),
  stringsAsFactors = FALSE
)
clinical <- multi[multi$Model == "Clinical_adjusted", ]
clinical$Model <- "Clinical adjusted"
clinical_caf <- multi[multi$Model == "Clinical_plus_CAF", ]
clinical_caf$Model <- "Clinical + CAF adjusted"

forest <- rbind(
  univ[, c("Cohort", "Model", "N", "Events", "HR", "CI_low", "CI_high", "P")],
  clinical[, c("Cohort", "Model", "N", "Events", "HR", "CI_low", "CI_high", "P")],
  clinical_caf[, c("Cohort", "Model", "N", "Events", "HR", "CI_low", "CI_high", "P")]
)
forest$I2_percent <- NA_real_
forest$Cohort <- sub("TCGA_PAAD", "TCGA-PAAD", forest$Cohort)
forest$Cohort <- sub("Moffitt_GSE71729", "Moffitt (GSE71729)", forest$Cohort)
forest$Cohort <- sub("Puleo_EMTAB6134", "Puleo (E-MTAB-6134)", forest$Cohort)

external_names <- c("Moffitt (GSE71729)", "Puleo (E-MTAB-6134)", "GSE57495")
meta_rows <- do.call(rbind, lapply(unique(forest$Model), function(model) {
  d <- forest[forest$Model == model & forest$Cohort %in% external_names, ]
  out <- fixed_meta(d)
  out$Model <- model
  out
}))
forest_all <- rbind(
  forest,
  meta_rows[, c(
    "Cohort", "Model", "N", "Events", "HR", "CI_low", "CI_high", "P",
    "I2_percent"
  )]
)
forest_all$Model <- factor(
  forest_all$Model,
  levels = c("Univariable", "Clinical adjusted", "Clinical + CAF adjusted")
)
forest_all$Cohort <- factor(
  forest_all$Cohort,
  levels = rev(c(
    "TCGA-PAAD", "Moffitt (GSE71729)", "Puleo (E-MTAB-6134)",
    "GSE57495", "External meta"
  ))
)
forest_all$Role <- ifelse(
  as.character(forest_all$Cohort) == "TCGA-PAAD", "Discovery",
  ifelse(as.character(forest_all$Cohort) == "External meta", "Meta", "Validation")
)
forest_all$Label <- sprintf(
  "%.2f (%.2f–%.2f), P=%s",
  forest_all$HR, forest_all$CI_low, forest_all$CI_high,
  format_p(forest_all$P)
)
is_meta <- forest_all$Role == "Meta"
forest_all$Label[is_meta] <- paste0(
  forest_all$Label[is_meta],
  sprintf("; I²=%.1f%%", forest_all$I2_percent[is_meta])
)
write.csv(
  forest_all,
  file.path(out_dir, "figure_forest_three_models_data.csv"),
  row.names = FALSE
)

role_colors <- c(
  Discovery = "#4D4D4D",
  Validation = "#0072B2",
  Meta = "#D55E00"
)
p_forest <- ggplot(
  forest_all,
  aes(x = HR, y = Cohort, color = Role, shape = Role)
) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55", linewidth = 0.45) +
  geom_segment(
    aes(x = CI_low, xend = CI_high, yend = Cohort),
    linewidth = 0.75
  ) +
  geom_point(size = 2.6, stroke = 0.8) +
  geom_text(
    aes(x = 3.35, label = Label),
    hjust = 0, size = 3.05, color = "black", show.legend = FALSE
  ) +
  facet_wrap(~Model, nrow = 1) +
  scale_x_log10(
    limits = c(0.75, 14),
    breaks = c(0.75, 1, 1.5, 2, 3, 5, 10)
  ) +
  scale_color_manual(values = role_colors) +
  scale_shape_manual(values = c(Discovery = 16, Validation = 16, Meta = 18)) +
  labs(
    x = "Hazard ratio per 1-SD increase in the 5-gene score (log scale)",
    y = NULL,
    color = NULL,
    shape = NULL,
    title = "Locked 5-gene remodeling–cytoskeleton score and overall survival",
    subtitle = "External meta combines Moffitt, Puleo and GSE57495; TCGA is excluded",
    caption = paste0(
      "Clinical adjustment: TCGA, age + AJCC stage + grade; Moffitt, T + N + margin + sex; ",
      "Puleo, T + N + grade + margin + sex; GSE57495, stage. ",
      "Clinical + CAF additionally includes the predefined CAF12 score."
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10.5) +
  theme(
    strip.background = element_rect(fill = "grey94", color = "grey70"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5),
    plot.caption = element_text(size = 8.2, hjust = 0),
    legend.position = "bottom",
    axis.text.y = element_text(size = 9),
    panel.spacing.x = unit(1.2, "lines"),
    plot.margin = margin(8, 15, 8, 8)
  )

ggsave(
  file.path(fig_dir, "Fig5A_four_cohort_adjusted_forest.png"),
  p_forest, width = 14.5, height = 6.2, units = "in", dpi = 300,
  bg = "white"
)
ggsave(
  file.path(fig_dir, "Fig5A_four_cohort_adjusted_forest.pdf"),
  p_forest, width = 14.5, height = 6.2, units = "in",
  device = cairo_pdf
)

# Leave-one-gene-out: cohort estimates plus external pooled estimates.
loo <- read.csv(
  file.path(out_dir, "five_gene_leave_one_out_by_cohort.csv"),
  stringsAsFactors = FALSE
)
loo <- loo[loo$Cohort %in% c(
  "Moffitt_GSE71729", "Puleo_EMTAB6134", "GSE57495"
), ]
loo$Cohort <- sub("Moffitt_GSE71729", "Moffitt", loo$Cohort)
loo$Cohort <- sub("Puleo_EMTAB6134", "Puleo", loo$Cohort)

loo_meta <- read.csv(
  file.path(out_dir, "five_gene_leave_one_out_external_meta.csv"),
  stringsAsFactors = FALSE
)
variant_labels <- c(
  Full5 = "Full 5-gene score",
  Drop_PLOD2 = "Without PLOD2",
  Drop_MSN = "Without MSN",
  Drop_TPM4 = "Without TPM4",
  Drop_PPFIBP1 = "Without PPFIBP1",
  Drop_SEC31A = "Without SEC31A"
)
variant_levels <- rev(unname(variant_labels))
loo$VariantLabel <- factor(variant_labels[loo$Variant], levels = variant_levels)
loo_meta$VariantLabel <- factor(
  variant_labels[loo_meta$Variant], levels = variant_levels
)

cohort_colors <- c(Moffitt = "#0072B2", Puleo = "#009E73", GSE57495 = "#CC79A7")
p_loo_cohort <- ggplot(
  loo,
  aes(x = HR, y = VariantLabel, color = Cohort)
) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55", linewidth = 0.45) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y", width = 0,
    position = position_dodge(width = 0.55), linewidth = 0.65
  ) +
  geom_point(position = position_dodge(width = 0.55), size = 2.2) +
  scale_x_log10(
    limits = c(0.75, 3.1),
    breaks = c(0.75, 1, 1.25, 1.5, 2, 3)
  ) +
  scale_color_manual(values = cohort_colors) +
  labs(
    x = "Cohort-specific HR per 1 SD (log scale)",
    y = NULL,
    color = NULL,
    title = "A  Individual validation cohorts"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "bottom"
  )

loo_meta$Label <- sprintf(
  "%.2f (%.2f–%.2f)\nP=%s; I²=%.1f%%",
  loo_meta$External_meta_HR, loo_meta$CI_low, loo_meta$CI_high,
  format_p(loo_meta$P), loo_meta$I2_percent
)
p_loo_meta <- ggplot(
  loo_meta,
  aes(x = External_meta_HR, y = VariantLabel)
) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55", linewidth = 0.45) +
  geom_segment(
    aes(x = CI_low, xend = CI_high, yend = VariantLabel),
    linewidth = 0.75, color = "#D55E00"
  ) +
  geom_point(size = 2.5, color = "#D55E00") +
  geom_text(
    aes(x = 1.72, label = Label), hjust = 0, size = 3.0
  ) +
  scale_x_log10(
    limits = c(0.95, 2.35),
    breaks = c(1, 1.25, 1.5, 2)
  ) +
  labs(
    x = "External fixed-effect meta HR (log scale)",
    y = NULL,
    title = "B  Leave-one-gene-out external meta-analysis"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(6, 90, 6, 6)
  )

p_loo <- p_loo_cohort + p_loo_meta +
  plot_layout(widths = c(1.08, 1)) +
  plot_annotation(
    title = "Stability of the locked 5-gene score",
    subtitle = "No leave-one-gene-out result was used to redefine the Full5 score",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9.5)
    )
  )

ggsave(
  file.path(fig_dir, "FigS_leave_one_gene_out_stability.png"),
  p_loo, width = 11.5, height = 5.8, units = "in", dpi = 300,
  bg = "white"
)
ggsave(
  file.path(fig_dir, "FigS_leave_one_gene_out_stability.pdf"),
  p_loo, width = 11.5, height = 5.8, units = "in",
  device = cairo_pdf
)

readme <- c(
  "# 冻结5基因评分统计图说明",
  paste0("Date: ", Sys.Date()),
  "",
  "## Fig5A_four_cohort_adjusted_forest",
  "- 效应量均为5基因评分每升高1 SD的OS HR。",
  "- Univariable：仅评分。",
  "- Clinical adjusted：TCGA=年龄+AJCC分期+分级；Moffitt=T+N+切缘+性别；Puleo=T+N+分级+切缘+性别；GSE57495=分期。",
  "- Clinical + CAF adjusted：在各自临床模型基础上加入预定义CAF12评分。",
  "- External meta仅合并Moffitt、Puleo、GSE57495，不包含TCGA。",
  "",
  "## FigS_leave_one_gene_out_stability",
  "- A：三个外部队列逐队列HR与95%CI。",
  "- B：三个外部队列固定效应meta及I²。",
  "- Full5为冻结主评分；Drop-one结果只用于判断是否由单一基因驱动。",
  "",
  "## 文件",
  "- PNG：300 dpi，便于预览和投稿排版。",
  "- PDF：矢量图，便于后续编辑。",
  "- figure_forest_three_models_data.csv：森林图完整底层数值。"
)
writeLines(readme, file.path(fig_dir, "FIGURE_README_CN.md"))

message("Saved locked-score figures to: ", fig_dir)
