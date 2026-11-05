suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))

root <- normalizePath(".", winslash = "/")
step18 <- file.path(root, "results", "step18_yellow317_mcode")
pkg <- file.path(
  root, "manuscript"
)
supp_dir <- file.path(pkg, "figures", "supplementary")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

locked <- read.csv(
  file.path(step18, "mcode1_four_blocks_locked.csv"),
  stringsAsFactors = FALSE
)
genes <- locked$Gene
edges <- read.csv(
  file.path(step18, "string_ppi_edges_score400.csv"),
  stringsAsFactors = FALSE
)
edges <- edges[
  edges$preferredName_A %in% genes & edges$preferredName_B %in% genes,
  c("preferredName_A", "preferredName_B", "score")
]
names(edges)[1:2] <- c("from", "to")
edges$key <- apply(edges[, c("from", "to")], 1, function(x) {
  paste(sort(x), collapse = "|")
})
edges <- do.call(rbind, lapply(split(edges, edges$key), function(d) {
  d[which.max(d$score), , drop = FALSE]
}))
g <- graph_from_data_frame(edges[, c("from", "to", "score")], directed = FALSE)
g <- simplify(
  g, remove.multiple = TRUE, remove.loops = TRUE,
  edge.attr.comb = list(score = "max")
)
E(g)$weight <- E(g)$score

truth <- setNames(as.integer(factor(locked$Block)), locked$Gene)
truth <- truth[V(g)$name]

set.seed(20260821)
methods <- list(
  Louvain = cluster_louvain(g, weights = E(g)$weight),
  Walktrap = cluster_walktrap(g, weights = E(g)$weight, steps = 4),
  FastGreedy = cluster_fast_greedy(g, weights = E(g)$weight),
  Leiden = cluster_leiden(
    g, objective_function = "modularity",
    weights = E(g)$weight, resolution = 1
  )
)

summary_rows <- do.call(rbind, lapply(names(methods), function(nm) {
  cl <- methods[[nm]]
  data.frame(
    Method = nm,
    Communities = length(unique(membership(cl))),
    Modularity = modularity(g, membership(cl), weights = E(g)$weight),
    ARI_vs_biology_blocks = compare(
      membership(cl), truth, method = "adjusted.rand"
    ),
    NMI_vs_biology_blocks = compare(
      membership(cl), truth, method = "nmi"
    )
  )
}))
print(summary_rows)

members <- data.frame(
  Gene = V(g)$name,
  BiologyBlock = locked$Block[match(V(g)$name, locked$Gene)],
  stringsAsFactors = FALSE
)
for (nm in names(methods)) {
  members[[nm]] <- membership(methods[[nm]])[members$Gene]
}
members <- members[order(members$BiologyBlock, members$Gene), ]
print(members, row.names = FALSE)

write.csv(
  summary_rows,
  file.path(step18, "four_block_unsupervised_community_agreement.csv"),
  row.names = FALSE
)
write.csv(
  members,
  file.path(step18, "four_block_unsupervised_membership.csv"),
  row.names = FALSE
)

# Visual sensitivity check: same PPI layout, biology annotation versus Walktrap.
set.seed(20260821)
xy <- layout_with_fr(g, weights = E(g)$weight, niter = 1500)
nodes <- data.frame(
  Gene = V(g)$name,
  x = xy[, 1],
  y = xy[, 2],
  BiologyBlock = locked$Block[match(V(g)$name, locked$Gene)],
  Walktrap = membership(methods$Walktrap)[V(g)$name],
  stringsAsFactors = FALSE
)

# Map each Walktrap community to its majority biology block only for readable
# labels and identification of boundary genes; clustering itself is unchanged.
community_map <- vapply(sort(unique(nodes$Walktrap)), function(k) {
  tab <- sort(
    table(nodes$BiologyBlock[nodes$Walktrap == k]),
    decreasing = TRUE
  )
  names(tab)[1]
}, character(1))
names(community_map) <- as.character(sort(unique(nodes$Walktrap)))
nodes$WalktrapLabel <- unname(community_map[as.character(nodes$Walktrap)])
nodes$Boundary <- nodes$BiologyBlock != nodes$WalktrapLabel

edge_plot <- as_data_frame(g, what = "edges")
edge_plot$x <- nodes$x[match(edge_plot$from, nodes$Gene)]
edge_plot$y <- nodes$y[match(edge_plot$from, nodes$Gene)]
edge_plot$xend <- nodes$x[match(edge_plot$to, nodes$Gene)]
edge_plot$yend <- nodes$y[match(edge_plot$to, nodes$Gene)]

block_labels <- c(
  Collagen_scaffold = "Collagen scaffold",
  Basement_ECM = "Basement membrane / ECM",
  Integrin_adhesion = "Integrin adhesion",
  Remodeling_cytoskeleton = "Remodeling / cytoskeleton"
)
block_colors <- c(
  Collagen_scaffold = "#0072B2",
  Basement_ECM = "#009E73",
  Integrin_adhesion = "#E69F00",
  Remodeling_cytoskeleton = "#CC79A7"
)
network_panel <- function(fill_var, title, subtitle) {
  ggplot() +
    geom_segment(
      data = edge_plot,
      aes(x = x, y = y, xend = xend, yend = yend, linewidth = score),
      color = "grey72", alpha = 0.65
    ) +
    geom_point(
      data = nodes[!nodes$Boundary, ],
      aes(x = x, y = y, fill = .data[[fill_var]]),
      shape = 21, size = 4.8, color = "white", stroke = 0.6
    ) +
    geom_point(
      data = nodes[nodes$Boundary, ],
      aes(x = x, y = y, fill = .data[[fill_var]]),
      shape = 21, size = 5.3, color = "black", stroke = 1.25
    ) +
    geom_text_repel(
      data = nodes,
      aes(x = x, y = y, label = Gene, fontface = ifelse(Boundary, "bold", "plain")),
      size = 2.8, max.overlaps = Inf, seed = 20260821,
      box.padding = 0.18, point.padding = 0.1, min.segment.length = 0
    ) +
    scale_fill_manual(
      values = block_colors, labels = block_labels,
      guide = guide_legend(nrow = 2)
    ) +
    scale_linewidth(range = c(0.2, 1.1), guide = "none") +
    coord_equal(clip = "off") +
    labs(title = title, subtitle = subtitle, fill = NULL) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8),
      legend.position = "bottom",
      plot.margin = margin(5, 10, 5, 10)
    )
}

p_biology <- network_panel(
  "BiologyBlock",
  "A  Biology-annotated blocks",
  "GO, known protein functions and network neighborhoods"
)
p_walktrap <- network_panel(
  "WalktrapLabel",
  "B  Unsupervised Walktrap communities",
  "4 communities; ARI=0.740, NMI=0.770; black outlines = boundary genes"
)

fig <- p_biology | p_walktrap
ggsave(
  file.path(supp_dir, "community_comparison_biology_walktrap.png"),
  fig, width = 12, height = 6.3, units = "in", dpi = 300, bg = "white"
)
ggsave(
  file.path(supp_dir, "community_comparison_biology_walktrap.pdf"),
  fig, width = 12, height = 6.3, units = "in", device = cairo_pdf
)
