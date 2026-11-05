#!/usr/bin/env Rscript

# Rebuild Figure 1 and Figure S1 directly from WGCNA plotting objects and
# export true vector PDF/SVG files (no embedded raster panels).

options(stringsAsFactors = FALSE)
set.seed(20260820)

suppressPackageStartupMessages({
  library(WGCNA)
  library(grid)
  library(gridGraphics)
  library(gridExtra)
  library(svglite)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_arg)) stop("Run with Rscript --file=...")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")

# Recreate the exact cleaned WGCNA objects used by the manuscript figures.
source(file.path(project_root, "code", "41_wgcna_figures_cleaned.R"), local = FALSE)

out_dir <- file.path(
  project_root, "manuscript",
  "figures", "results_section1_ordered"
)
main_dir <- file.path(
  project_root, "manuscript",
  "figures", "main"
)
supp_dir <- file.path(
  project_root, "manuscript",
  "figures", "supplementary"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set_text_size <- function(grob, size = 8) {
  grob <- grid.force(grob)
  if (inherits(grob, "text")) {
    if (is.null(grob$gp)) grob$gp <- gpar()
    grob$gp$fontsize <- size
  }
  if (!is.null(grob$children) && length(grob$children)) {
    grob$children <- do.call(
      gList,
      lapply(as.list(grob$children), set_text_size, size = size)
    )
  }
  if (!is.null(grob$grobs) && length(grob$grobs)) {
    grob$grobs <- lapply(grob$grobs, set_text_size, size = size)
  }
  grob
}

capture_base_grob <- function(draw, width, height) {
  pdf(NULL, width = width, height = height, pointsize = 8)
  on.exit(dev.off(), add = TRUE)
  # File devices keep the display list off by default; grid.echo() needs it
  # enabled to replay base/WGCNA graphics as editable vector grobs.
  dev.control(displaylist = "enable")
  draw()
  grid.echo()
  set_text_size(grid.grab(), size = 8)
}

panel_label <- function(grob, label) {
  grobTree(
    grob,
    textGrob(
      label, x = unit(0.012, "npc"), y = unit(0.988, "npc"),
      just = c("left", "top"),
      gp = gpar(fontface = "bold", fontsize = 10)
    )
  )
}

draw_soft_threshold <- function() {
  par(
    mfrow = c(1, 2), mar = c(5, 5, 3.5, 1),
    ps = 8, cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1
  )
  plot(
    fit$Power, fit$SFT.R.sq,
    xlab = "Soft threshold (power)",
    ylab = "Scale-free fit, signed R2",
    type = "n", main = "Scale independence"
  )
  text(fit$Power, fit$SFT.R.sq, labels = fit$Power, cex = 1, col = "firebrick")
  abline(h = 0.80, col = "dodgerblue", lty = 2, lwd = 1.5)
  points(
    power, fit$SFT.R.sq[fit$Power == power],
    pch = 21, bg = "forestgreen", cex = 1.8
  )
  plot(
    fit$Power, fit$mean.k.,
    xlab = "Soft threshold (power)", ylab = "Mean connectivity",
    type = "n", main = "Mean connectivity"
  )
  text(fit$Power, fit$mean.k., labels = fit$Power, cex = 1, col = "firebrick")
  points(
    power, fit$mean.k.[fit$Power == power],
    pch = 21, bg = "forestgreen", cex = 1.8
  )
}

draw_gene_tree <- function() {
  # plotDendroAndColors() drives its own layout(), and the axis titles it passes
  # through do not survive grid.echo(); draw them in the outer margin instead.
  par(
    mar = c(5, 4, 3, 2), oma = c(1.2, 2.2, 0, 0),
    ps = 8, cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1
  )
  plotDendroAndColors(
    geneTree, moduleColors[net$blockGenes[[1]]], "Module colors",
    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
    main = "Gene dendrogram and modules",
    xlab = "",
    ylab = "",
    sub = "",
    cex.colorLabels = 1
  )
  mtext("Height", side = 2, line = 0.6, outer = TRUE, cex = 1)
  mtext("Genes", side = 1, line = 0.1, outer = TRUE, cex = 1)
}

draw_module_trait <- function() {
  par(
    mar = c(7, 11, 3, 2),
    ps = 8, cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1
  )
  trait_labels <- c(
    CAF_score = "CAF score",
    OS_time = "OS time",
    OS_event = "OS event",
    Age = "Age",
    Stage = "Stage",
    Purity = "Purity"
  )
  x_labs <- unname(trait_labels[colnames(traits)])
  x_labs[is.na(x_labs)] <- colnames(traits)[is.na(x_labs)]
  labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = x_labs,
    yLabels = rownames(moduleTraitCor),
    ySymbols = rownames(moduleTraitCor),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 1,
    cex.lab = 1,
    zlim = c(-1, 1),
    main = "Module-trait correlations"
  )
}

draw_mm_gs <- function() {
  par(
    mar = c(5, 5, 3, 2),
    ps = 8, cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1
  )
  # Hub genes crowd the top-right corner, so labels are moved into a right-hand
  # gutter and spread vertically with leader lines instead of overplotting.
  x_range <- range(mm)
  gutter <- 0.22 * diff(x_range)
  plot(
    mm, gs,
    xlab = paste0("Module membership |MM| in ", best_color),
    ylab = "Gene significance |GS| for CAF score",
    main = "Hub gene identification",
    xlim = c(x_range[1], x_range[2] + gutter),
    pch = 21, bg = "darkgreen", col = "darkgreen", cex = 0.85
  )
  abline(v = 0.8, lty = 2, col = "grey40", lwd = 1.2)
  abline(h = 0.3, lty = 2, col = "grey40", lwd = 1.2)
  hub_idx <- which(mm > 0.8 & gs > 0.3)
  if (length(hub_idx)) {
    ord <- order(-(mm[hub_idx] * gs[hub_idx]))
    top <- head(hub_idx[ord], 8)
    top <- top[order(-gs[top])]
    label_x <- x_range[2] + 0.05 * diff(x_range)
    label_y <- gs[top]
    min_gap <- 0.042 * diff(range(gs))
    for (i in seq_along(label_y)[-1]) {
      if (label_y[i - 1] - label_y[i] < min_gap) {
        label_y[i] <- label_y[i - 1] - min_gap
      }
    }
    segments(
      mm[top], gs[top], label_x - 0.004, label_y,
      col = "grey60", lwd = 0.6
    )
    text(label_x, label_y, labels = gnames[top], cex = 1, pos = 4, offset = 0)
  }
}

draw_before_qc <- function() {
  par(
    mar = c(5, 4, 3, 2),
    ps = 8, cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1
  )
  trait_labels <- c(
    CAF_score = "CAF score",
    OS_time = "OS time",
    OS_event = "OS event",
    Age = "Age",
    Stage = "Stage",
    Purity = "Purity"
  )
  group_labs0 <- unname(trait_labels[colnames(traits0)])
  group_labs0[is.na(group_labs0)] <- colnames(traits0)[is.na(group_labs0)]
  plotDendroAndColors(
    sampleTree0,
    numbers2colors(traits0, signed = FALSE, naColor = "grey80"),
    groupLabels = group_labs0,
    main = paste0("Before outlier removal (n=", nrow(datExpr0), ")"),
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05, cex.colorLabels = 1,
    xlab = "Samples",
    ylab = "Height",
    sub = ""
  )
  abline(h = cut_h, col = "red", lty = 2, lwd = 1.5)
}

draw_after_qc <- function() {
  par(
    mar = c(5, 4, 3, 2),
    ps = 8, cex = 1, cex.axis = 1, cex.lab = 1, cex.main = 1
  )
  trait_labels <- c(
    CAF_score = "CAF score",
    OS_time = "OS time",
    OS_event = "OS event",
    Age = "Age",
    Stage = "Stage",
    Purity = "Purity"
  )
  group_labs <- unname(trait_labels[colnames(traits)])
  group_labs[is.na(group_labs)] <- colnames(traits)[is.na(group_labs)]
  plotDendroAndColors(
    sampleTree,
    numbers2colors(traits, signed = FALSE, naColor = "grey80"),
    groupLabels = group_labs,
    main = paste0("After outlier removal (n=", nrow(datExpr), ")"),
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05, cex.colorLabels = 1,
    xlab = "Samples",
    ylab = "Height",
    sub = ""
  )
}

# Fig.1C: show Pearson r for all cells; annotate significance with stars
# (* P<0.05, ** P<0.01, *** P<0.001) instead of embedding raw P values.
textMatrix <- matrix("", nrow = nrow(moduleTraitCor), ncol = ncol(moduleTraitCor))
for (i in seq_len(nrow(moduleTraitCor))) {
  for (j in seq_len(ncol(moduleTraitCor))) {
    r_lab <- formatC(moduleTraitCor[i, j], digits = 2, format = "f")
    p_ij <- moduleTraitPvalue[i, j]
    stars <- if (is.na(p_ij)) {
      ""
    } else if (p_ij < 0.001) {
      "***"
    } else if (p_ij < 0.01) {
      "**"
    } else if (p_ij < 0.05) {
      "*"
    } else {
      ""
    }
    textMatrix[i, j] <- if (nzchar(stars)) paste0(r_lab, "\n", stars) else r_lab
  }
}

g_a <- panel_label(capture_base_grob(draw_soft_threshold, 10, 5), "A")
g_b <- panel_label(capture_base_grob(draw_gene_tree, 11, 7), "B")
g_c <- panel_label(capture_base_grob(draw_module_trait, 9.5, 10), "C")
g_d <- panel_label(capture_base_grob(draw_mm_gs, 7.5, 7), "D")
fig1 <- arrangeGrob(
  g_a, g_b, g_c, g_d,
  ncol = 2,
  widths = c(1, 1),
  heights = c(0.82, 1.18)
)

g_s1a <- panel_label(capture_base_grob(draw_before_qc, 11, 7), "A")
g_s1b <- panel_label(capture_base_grob(draw_after_qc, 11, 7), "B")
fig_s1 <- arrangeGrob(
  g_s1a, g_s1b, ncol = 1,
  bottom = textGrob(
    paste(
      "Fig. S1 | Sample quality control for WGCNA discovery cohort.",
      paste(
        "Hierarchical clustering of TCGA-PAAD samples based on expression",
        "dissimilarity, with trait annotations (CAF score, OS time, OS event)."
      ),
      "(A) Before outlier removal (n=133).",
      paste(
        "(B) After removing outlier TCGA-2J-AABV-01A, yielding the final",
        "analytical cohort (n=132)."
      ),
      sep = "\n"
    ),
    x = 0.01, hjust = 0, vjust = 0.5,
    gp = gpar(fontsize = 8, lineheight = 1.15)
  )
)

write_vector_pair <- function(grob, stem, width, height) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(out_dir, paste0(stem, ".svg"))
  tmp_pdf <- tempfile(fileext = ".pdf")
  tmp_svg <- tempfile(fileext = ".svg")

  ok_pdf <- FALSE
  tryCatch({
    cairo_pdf(tmp_pdf, width = width, height = height, family = "sans", pointsize = 8)
    grid.newpage()
    grid.draw(grob)
    dev.off()
    ok_pdf <- TRUE
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
  })
  if (!ok_pdf) {
    pdf(tmp_pdf, width = width, height = height, useDingbats = FALSE, pointsize = 8)
    grid.newpage()
    grid.draw(grob)
    dev.off()
  }

  svglite(
    tmp_svg, width = width, height = height,
    system_fonts = list(sans = "Arial"), pointsize = 8
  )
  grid.newpage()
  grid.draw(grob)
  dev.off()

  # svglite uses editable <text> elements, but nested WGCNA/base viewports can
  # rescale some fonts. Normalize the delivered SVG explicitly to journal
  # typography: all plot text 8 pt and panel letters 10 pt. Remove textLength
  # so editors retain natural editable Arial text rather than stretched glyphs.
  svg_lines <- readLines(tmp_svg, warn = FALSE, encoding = "UTF-8")
  is_text <- grepl("<text\\b", svg_lines)
  is_panel_letter <- is_text &
    grepl(">[A-G]</text>", svg_lines) &
    grepl("font-weight: bold|font-weight: 700|font-weight:700", svg_lines)
  svg_lines[is_text & !is_panel_letter] <- gsub(
    "font-size: [0-9.]+px", "font-size: 8.00pt",
    svg_lines[is_text & !is_panel_letter]
  )
  svg_lines[is_panel_letter] <- gsub(
    "font-size: [0-9.]+px", "font-size: 10.00pt",
    svg_lines[is_panel_letter]
  )
  svg_lines[is_text] <- gsub(
    " textLength='[^']*' lengthAdjust='spacingAndGlyphs'",
    "", svg_lines[is_text]
  )
  writeLines(svg_lines, tmp_svg, useBytes = TRUE)

  if (!file.exists(tmp_pdf) || file.info(tmp_pdf)$size < 10000) {
    stop("Invalid PDF export: ", tmp_pdf)
  }
  if (!file.exists(tmp_svg) || file.info(tmp_svg)$size < 10000) {
    stop("Invalid SVG export: ", tmp_svg)
  }
  file.copy(tmp_pdf, pdf_path, overwrite = TRUE)
  file.copy(tmp_svg, svg_path, overwrite = TRUE)
  unlink(c(tmp_pdf, tmp_svg))
}

stem1 <- "01_Figure1_CAF_anchored_yellow_module"
stems1 <- "02_FigureS1_WGCNA_sample_QC"
write_vector_pair(fig1, stem1, width = 14, height = 11)
write_vector_pair(fig_s1, stems1, width = 12, height = 14.4)

file.copy(
  file.path(out_dir, paste0(stem1, ".pdf")),
  file.path(main_dir, "Figure1_lock_yellow_CAF_module.pdf"),
  overwrite = TRUE
)
file.copy(
  file.path(out_dir, paste0(stem1, ".svg")),
  file.path(main_dir, "Figure1_lock_yellow_CAF_module.svg"),
  overwrite = TRUE
)
file.copy(
  file.path(out_dir, paste0(stems1, ".pdf")),
  file.path(supp_dir, "FigureS1_WGCNA_construction_QC.pdf"),
  overwrite = TRUE
)
file.copy(
  file.path(out_dir, paste0(stems1, ".svg")),
  file.path(supp_dir, "FigureS1_WGCNA_construction_QC.svg"),
  overwrite = TRUE
)
png_s1 <- file.path(out_dir, paste0(stems1, ".png"))
png(png_s1, width = 12, height = 14.4, units = "in", res = 300, bg = "white")
grid.newpage()
grid.draw(fig_s1)
dev.off()
file.copy(png_s1, file.path(supp_dir, "FigureS1_WGCNA_construction_QC.png"), overwrite = TRUE)
sub_supp <- file.path(
  project_root, "manuscript",
  "submission_pack", "figures", "supplementary"
)
dir.create(sub_supp, recursive = TRUE, showWarnings = FALSE)
for (ext in c("pdf", "svg", "png")) {
  file.copy(
    file.path(supp_dir, paste0("FigureS1_WGCNA_construction_QC.", ext)),
    file.path(sub_supp, paste0("FigureS1_WGCNA_construction_QC.", ext)),
    overwrite = TRUE
  )
}

message("True-vector PDF/SVG exports written to: ", out_dir)
