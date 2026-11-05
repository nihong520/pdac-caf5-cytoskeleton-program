#!/usr/bin/env Rscript
# Export Fig.2, Fig.3, Fig.S2, Fig.S3 (+ S6G panel) as true-vector PDF/SVG
# matching Fig.1 typography: editable SVG text, 8 pt body, 10 pt panel tags.

options(stringsAsFactors = FALSE)
set.seed(20260821)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(igraph)
  library(ggrepel)
  library(svglite)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(script_arg)) stop("Run with Rscript --file=...")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")

pkg <- file.path(root, "manuscript")
main_dir <- file.path(pkg, "figures", "main")
supp_dir <- file.path(pkg, "figures", "supplementary")
out_dir <- file.path(pkg, "figures", "results_section2_4_ordered")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

step18 <- file.path(root, "results", "step18_yellow317_mcode")
step19 <- file.path(root, "results", "step19_five_gene_locked_validation")
step20 <- file.path(root, "results", "step20_statistical_strengthening")
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

theme_journal_void <- function(base = 8) {
  theme_void(base_size = base) %+replace%
    theme(
      plot.title = element_text(face = "bold", size = 8, hjust = 0),
      plot.subtitle = element_text(size = 8, hjust = 0),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      plot.tag = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      complete = TRUE
    )
}

normalize_svg_text <- function(svg_path) {
  svg_lines <- readLines(svg_path, warn = FALSE, encoding = "UTF-8")
  is_text <- grepl("<text\\b", svg_lines)
  # Panel tags from patchwork are often bold single letters near the start of
  # each panel; keep explicit 10 pt for short bold A-G tags.
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
    file.copy(pdf_path, file.path(main_dir, paste0(also_main, ".pdf")), overwrite = TRUE)
    file.copy(svg_path, file.path(main_dir, paste0(also_main, ".svg")), overwrite = TRUE)
    file.copy(png_path, file.path(main_dir, paste0(also_main, ".png")), overwrite = TRUE)
  }
  if (!is.null(also_supp)) {
    file.copy(pdf_path, file.path(supp_dir, paste0(also_supp, ".pdf")), overwrite = TRUE)
    file.copy(svg_path, file.path(supp_dir, paste0(also_supp, ".svg")), overwrite = TRUE)
    file.copy(png_path, file.path(supp_dir, paste0(also_supp, ".png")), overwrite = TRUE)
    sub_supp <- file.path(pkg, "submission_pack", "figures", "supplementary")
    dir.create(sub_supp, recursive = TRUE, showWarnings = FALSE)
    file.copy(pdf_path, file.path(sub_supp, paste0(also_supp, ".pdf")), overwrite = TRUE)
    file.copy(svg_path, file.path(sub_supp, paste0(also_supp, ".svg")), overwrite = TRUE)
    file.copy(png_path, file.path(sub_supp, paste0(also_supp, ".png")), overwrite = TRUE)
  }
  message("Wrote ", stem)
}

# ---------------------------------------------------------------------------
# Figure 2 panels (content aligned with manuscript: A PPI, B GO, C blocks, D meta)
# ---------------------------------------------------------------------------
cluster_summary <- read.csv(
  file.path(step18, "mcode_clusters_summary.csv"), stringsAsFactors = FALSE
)
mcode1_genes <- strsplit(
  cluster_summary$Genes[cluster_summary$Cluster == "MCODE1"], ";"
)[[1]]
edges <- read.csv(
  file.path(step18, "string_ppi_edges_score400.csv"), stringsAsFactors = FALSE
)
edges <- edges[
  edges$preferredName_A %in% mcode1_genes &
    edges$preferredName_B %in% mcode1_genes,
]
edge_df <- data.frame(
  from = edges$preferredName_A,
  to = edges$preferredName_B,
  score = edges$score
)
edge_df$key <- apply(edge_df[, c("from", "to")], 1, function(x) {
  paste(sort(x), collapse = "|")
})
edge_df <- do.call(rbind, lapply(split(edge_df, edge_df$key), function(d) {
  d[which.max(d$score), , drop = FALSE]
}))
g <- graph_from_data_frame(edge_df[, c("from", "to", "score")], directed = FALSE)
xy <- layout_with_fr(g, weights = E(g)$score, niter = 1500)
nodes <- data.frame(
  Gene = V(g)$name, x = xy[, 1], y = xy[, 2], Degree = degree(g)
)
categories <- list(
  "Collagen scaffold" = c(
    "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A2",
    "COL6A3", "COL8A1", "COL10A1", "COL11A1"
  ),
  "Basement membrane / ECM" = c(
    "FN1", "HSPG2", "LAMA4", "LAMB1", "LAMC1", "NID1", "NID2"
  ),
  "Integrin adhesion" = c("ITGAV", "ITGB1", "ITGB5", "PARVA", "VCL"),
  "Cytoskeleton remodeling" = c("PLOD2", "MSN", "TPM4", "PPFIBP1", "SEC31A")
)
nodes$Category <- NA_character_
for (nm in names(categories)) {
  nodes$Category[nodes$Gene %in% categories[[nm]]] <- nm
}
edge_coords <- as_data_frame(g, what = "edges")
edge_coords$x <- nodes$x[match(edge_coords$from, nodes$Gene)]
edge_coords$y <- nodes$y[match(edge_coords$from, nodes$Gene)]
edge_coords$xend <- nodes$x[match(edge_coords$to, nodes$Gene)]
edge_coords$yend <- nodes$y[match(edge_coords$to, nodes$Gene)]
category_colors <- c(
  "Collagen scaffold" = "#0072B2",
  "Basement membrane / ECM" = "#009E73",
  "Integrin adhesion" = "#E69F00",
  "Cytoskeleton remodeling" = "#CC79A7"
)

p_network <- ggplot() +
  geom_segment(
    data = edge_coords,
    aes(x = x, y = y, xend = xend, yend = yend, linewidth = score),
    color = "grey70", alpha = 0.6
  ) +
  geom_point(
    data = nodes,
    aes(x = x, y = y, size = Degree),
    shape = 21, fill = "#6BAED6", color = "white", stroke = 0.5
  ) +
  geom_text_repel(
    data = nodes,
    aes(x = x, y = y, label = Gene),
    size = 2.5, max.overlaps = Inf, seed = 20260821,
    box.padding = 0.2, point.padding = 0.12, min.segment.length = 0
  ) +
  scale_linewidth(range = c(0.2, 1.1), guide = "none") +
  scale_size_continuous(range = c(3, 7.5), guide = "none") +
  coord_equal(clip = "off") +
  labs(
    title = "Dense 26-gene PPI core"
  ) +
  theme_journal_void()

go <- read.csv(file.path(step18, "MCODE1_GO_BP.csv"), stringsAsFactors = FALSE)
keep_terms <- c(
  "extracellular matrix organization",
  "basement membrane organization",
  "collagen fibril organization",
  "integrin-mediated signaling pathway",
  "cell-substrate adhesion",
  "cell-matrix adhesion"
)
go_sel <- go[match(keep_terms, go$Description), ]
go_sel <- go_sel[complete.cases(go_sel$Description), ]
go_sel$Term <- vapply(
  go_sel$Description,
  function(x) paste(strwrap(x, width = 32), collapse = "\n"),
  character(1)
)
go_sel$Term <- factor(go_sel$Term, levels = rev(go_sel$Term))
go_sel$MinusLog10FDR <- -log10(go_sel$p.adjust)

p_go <- ggplot(go_sel, aes(x = RichFactor, y = Term)) +
  geom_segment(
    aes(x = 0, xend = RichFactor, yend = Term),
    linewidth = 1.2, color = "grey82", lineend = "round"
  ) +
  geom_point(
    aes(size = Count, color = MinusLog10FDR),
    alpha = 0.98
  ) +
  geom_text(
    aes(label = paste0("n=", Count)),
    nudge_x = 0.014, hjust = 0, size = 2.3, color = "grey25"
  ) +
  scale_x_continuous(
    limits = c(0, max(go_sel$RichFactor) + 0.055),
    breaks = seq(0, 0.30, by = 0.10),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_viridis_c(option = "C", direction = -1) +
  scale_size_continuous(range = c(4, 8)) +
  labs(
    title = "Functional enrichment",
    x = "Rich factor", y = NULL,
    color = expression(-log[10]("FDR")), size = "n genes"
  ) +
  theme_journal() +
  theme(
    axis.text.y = element_text(size = 8, lineheight = 0.95),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

block_cards <- data.frame(
  Category = names(categories),
  Title = c(
    "Collagen scaffold  (n=9)",
    "Basement membrane / ECM  (n=7)",
    "Integrin adhesion  (n=5)",
    "Cytoskeleton remodeling  (n=5)"
  ),
  Function = c(
    "Fibrillar and network-forming collagen structure",
    "Basement-membrane assembly and matrix organization",
    "Cell–matrix sensing and focal-adhesion coupling",
    "Matrix maturation, cytoskeletal remodeling and secretion"
  ),
  Genes = vapply(categories, function(x) {
    paste(strwrap(paste(x, collapse = ", "), width = 38), collapse = "\n")
  }, character(1)),
  xmin = c(0.00, 1.15, 0.00, 1.15),
  xmax = c(1.05, 2.20, 1.05, 2.20),
  ymin = c(1.15, 1.15, 0.00, 0.00),
  ymax = c(2.15, 2.15, 1.00, 1.00),
  stringsAsFactors = FALSE
)
block_cards$x <- block_cards$xmin + 0.07
block_cards$title_y <- block_cards$ymax - 0.13
block_cards$function_y <- block_cards$ymax - 0.37
block_cards$genes_y <- block_cards$ymax - 0.61
p_arch <- ggplot(block_cards) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Category),
    alpha = 0.14, color = NA
  ) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = Category),
    fill = NA, linewidth = 1
  ) +
  geom_text(
    aes(x = x, y = title_y, label = Title, color = Category),
    hjust = 0, vjust = 1, size = 2.8, fontface = "bold"
  ) +
  geom_text(
    aes(x = x, y = function_y, label = Function),
    hjust = 0, vjust = 1, size = 2.4, color = "grey25"
  ) +
  geom_text(
    aes(x = x, y = genes_y, label = Genes),
    hjust = 0, vjust = 1, size = 2.4, lineheight = 1.05
  ) +
  scale_fill_manual(values = category_colors, guide = "none") +
  scale_color_manual(values = category_colors, guide = "none") +
  coord_cartesian(xlim = c(-0.03, 2.23), ylim = c(-0.03, 2.22), clip = "off") +
  labs(
    title = "Four functional blocks"
  ) +
  theme_journal_void()

block_meta <- read.csv(
  file.path(step18, "mcode1_four_blocks_external_meta.csv"),
  stringsAsFactors = FALSE
)
block_labels <- c(
  Collagen_scaffold = "Collagen scaffold",
  Basement_ECM = "Basement membrane / ECM",
  Integrin_adhesion = "Integrin adhesion",
  Remodeling_cytoskeleton = "Cytoskeleton remodeling (locked)"
)
block_meta$BlockLabel <- factor(
  block_labels[block_meta$Block], levels = rev(unname(block_labels))
)
block_meta$Model <- factor(
  block_meta$Model,
  levels = c("Univariable", "CAF_adjusted"),
  labels = c("Univariable", "CAF adjusted")
)
p_blocks <- ggplot(
  block_meta,
  aes(x = External_meta_HR, y = BlockLabel, color = Model, shape = Model)
) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y", width = 0,
    position = position_dodge(width = 0.55), linewidth = 0.7
  ) +
  geom_point(position = position_dodge(width = 0.55), size = 2.3) +
  scale_x_log10(limits = c(0.7, 2.4), breaks = c(0.75, 1, 1.5, 2)) +
  scale_color_manual(values = c("Univariable" = "#0072B2", "CAF adjusted" = "#D55E00")) +
  scale_shape_manual(values = c("Univariable" = 16, "CAF adjusted" = 17)) +
  labs(
    title = "External comparison and locking",
    x = "Meta HR per 1 SD (log)", y = NULL, color = NULL, shape = NULL
  ) +
  theme_journal() +
  theme(
    axis.text.y = element_text(
      face = ifelse(
        levels(block_meta$BlockLabel) == "Cytoskeleton remodeling (locked)",
        "bold", "plain"
      )
    )
  )

fig2 <- (p_network | p_go) / (p_arch | p_blocks) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 10, face = "bold"))
  ) &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_vector(
  fig2, "01_Figure2_functional_dissection_locking",
  width = 14.5, height = 11,
  also_main = "Figure2_functional_dissection_locking"
)

# ---------------------------------------------------------------------------
# Figure S2: biology vs Walktrap
# ---------------------------------------------------------------------------
locked <- read.csv(
  file.path(step18, "mcode1_four_blocks_locked.csv"), stringsAsFactors = FALSE
)
genes <- locked$Gene
edges2 <- read.csv(
  file.path(step18, "string_ppi_edges_score400.csv"), stringsAsFactors = FALSE
)
edges2 <- edges2[
  edges2$preferredName_A %in% genes & edges2$preferredName_B %in% genes,
  c("preferredName_A", "preferredName_B", "score")
]
names(edges2)[1:2] <- c("from", "to")
edges2$key <- apply(edges2[, c("from", "to")], 1, function(x) paste(sort(x), collapse = "|"))
edges2 <- do.call(rbind, lapply(split(edges2, edges2$key), function(d) d[which.max(d$score), , drop = FALSE]))
g2 <- graph_from_data_frame(edges2[, c("from", "to", "score")], directed = FALSE)
g2 <- simplify(g2, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(score = "max"))
E(g2)$weight <- E(g2)$score
wt <- cluster_walktrap(g2, weights = E(g2)$weight, steps = 4)
xy2 <- layout_with_fr(g2, weights = E(g2)$weight, niter = 1500)
nodes2 <- data.frame(
  Gene = V(g2)$name, x = xy2[, 1], y = xy2[, 2],
  BiologyBlock = locked$Block[match(V(g2)$name, locked$Gene)],
  Walktrap = membership(wt)[V(g2)$name],
  stringsAsFactors = FALSE
)
community_map <- vapply(sort(unique(nodes2$Walktrap)), function(k) {
  tab <- sort(table(nodes2$BiologyBlock[nodes2$Walktrap == k]), decreasing = TRUE)
  names(tab)[1]
}, character(1))
names(community_map) <- as.character(sort(unique(nodes2$Walktrap)))
nodes2$WalktrapLabel <- unname(community_map[as.character(nodes2$Walktrap)])
nodes2$Boundary <- nodes2$BiologyBlock != nodes2$WalktrapLabel
nodes2$WalktrapCommunity <- factor(
  paste0("Community ", nodes2$Walktrap),
  levels = paste0("Community ", sort(unique(nodes2$Walktrap)))
)
edge_plot <- as_data_frame(g2, what = "edges")
edge_plot$x <- nodes2$x[match(edge_plot$from, nodes2$Gene)]
edge_plot$y <- nodes2$y[match(edge_plot$from, nodes2$Gene)]
edge_plot$xend <- nodes2$x[match(edge_plot$to, nodes2$Gene)]
edge_plot$yend <- nodes2$y[match(edge_plot$to, nodes2$Gene)]
block_colors_s2 <- c(
  Collagen_scaffold = "#0072B2",
  Basement_ECM = "#009E73",
  Integrin_adhesion = "#E69F00",
  Remodeling_cytoskeleton = "#CC79A7"
)
block_labels_s2 <- c(
  Collagen_scaffold = "Collagen scaffold",
  Basement_ECM = "Basement membrane / ECM",
  Integrin_adhesion = "Integrin adhesion",
  Remodeling_cytoskeleton = "Cytoskeleton remodeling"
)
walktrap_colors <- setNames(
  unname(block_colors_s2[community_map[as.character(sort(unique(nodes2$Walktrap)))]]),
  paste0("Community ", sort(unique(nodes2$Walktrap)))
)
network_panel <- function(fill_var, title, fill_values, fill_labels = NULL) {
  p <- ggplot() +
    geom_segment(
      data = edge_plot,
      aes(x = x, y = y, xend = xend, yend = yend, linewidth = score),
      color = "grey72", alpha = 0.65
    ) +
    geom_point(
      data = nodes2[!nodes2$Boundary, ],
      aes(x = x, y = y, fill = .data[[fill_var]]),
      shape = 21, size = 4.5, color = "white", stroke = 0.6
    ) +
    geom_point(
      data = nodes2[nodes2$Boundary, ],
      aes(x = x, y = y, fill = .data[[fill_var]]),
      shape = 21, size = 5.0, color = "black", stroke = 1.2
    ) +
    geom_text_repel(
      data = nodes2,
      aes(x = x, y = y, label = Gene, fontface = ifelse(Boundary, "bold", "plain")),
      size = 2.5, max.overlaps = Inf, seed = 20260821,
      box.padding = 0.18, point.padding = 0.1, min.segment.length = 0
    ) +
    scale_linewidth(range = c(0.2, 1.1), guide = "none") +
    coord_equal(clip = "off") +
    labs(title = title, fill = NULL) +
    theme_journal_void()
  if (is.null(fill_labels)) {
    p <- p + scale_fill_manual(values = fill_values, guide = guide_legend(nrow = 2))
  } else {
    p <- p + scale_fill_manual(
      values = fill_values, labels = fill_labels, guide = guide_legend(nrow = 2)
    )
  }
  p
}
p_s2a <- network_panel(
  "BiologyBlock",
  "Biology-annotated blocks",
  block_colors_s2,
  block_labels_s2
)
p_s2b <- network_panel(
  "WalktrapCommunity",
  "Walktrap communities",
  walktrap_colors
)
fig_s2 <- (p_s2a | p_s2b) +
  plot_annotation(
    tag_levels = "A",
    caption = paste(
      "Fig. S2 | Topology sensitivity of biology-annotated blocks versus Walktrap communities.",
      "The same 26-gene PPI network colored by biology-annotated blocks (A) or Walktrap communities (B).",
      "Walktrap yielded four communities (ARI=0.740, NMI=0.770); black outlines mark boundary genes (PLOD2, FN1, VCL).",
      sep = "\n"
    ),
    theme = theme(
      plot.tag = element_text(size = 10, face = "bold"),
      plot.caption = element_text(size = 8, hjust = 0, lineheight = 1.2, margin = margin(t = 8))
    )
  ) &
  theme(plot.tag = element_text(size = 10, face = "bold"))
save_vector(
  fig_s2, "02_FigureS2_PPI_community_sensitivity",
  width = 12, height = 7.6,
  also_supp = "FigureS2_PPI_community_sensitivity"
)

# Fig. 3 / S3 are owned by code/75_export_fig3_with_panel_f.R — do not overwrite.
