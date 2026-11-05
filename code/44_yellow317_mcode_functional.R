#!/usr/bin/env Rscript
# Yellow317: STRING PPI -> MCODE-like dense clusters -> per-cluster GO/KEGG
# Framing: functional dissection of CAF-ECM program (Option 1 supplement).
# Prognosis on top clusters is ONE locked side-check only (not signature fishing).

suppressPackageStartupMessages({
  library(igraph)
  library(httr)
  library(jsonlite)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

project_root <- normalizePath(getwd())
out_dir <- file.path(project_root, "results", "step18_yellow317_mcode")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gene_file <- file.path(
  project_root, "results", "step14_caf_wgcna",
  "wgcna_standard_figures", "selected_module_genes_AFTER_QC.csv"
)
genes <- unique(toupper(read.csv(gene_file, stringsAsFactors = FALSE)$Gene))
message("Yellow genes: ", length(genes))

# ---- STRING PPI (REST) ----
string_post <- function(path, body_list) {
  url <- paste0("https://string-db.org/api/", path)
  r <- httr::POST(url, body = body_list, encode = "form", httr::timeout(180))
  if (httr::status_code(r) >= 400) {
    stop("STRING API failed: ", httr::status_code(r), " ", httr::content(r, "text"))
  }
  httr::content(r, as = "text", encoding = "UTF-8")
}

ids_txt <- string_post("json/get_string_ids", list(
  identifiers = paste(genes, collapse = "\r"),
  species = 9606,
  limit = 1,
  echo_query = 1
))
id_map <- jsonlite::fromJSON(ids_txt)
if (is.null(nrow(id_map)) && length(id_map)) id_map <- as.data.frame(id_map)
id_map$queryItem <- toupper(as.character(id_map$queryItem))
id_map$stringId <- as.character(id_map$stringId)
id_map$preferredName <- toupper(as.character(id_map$preferredName))
write.csv(id_map, file.path(out_dir, "string_id_map.csv"), row.names = FALSE)
message("Mapped to STRING: ", nrow(id_map), " / ", length(genes))

net_txt <- string_post("tsv/network", list(
  identifiers = paste(id_map$stringId, collapse = "\r"),
  species = 9606,
  required_score = 400,
  network_type = "physical"
))
net <- utils::read.delim(textConnection(net_txt), stringsAsFactors = FALSE)
if (!nrow(net)) stop("No STRING edges at score>=400")
net$preferredName_A <- toupper(net$preferredName_A)
net$preferredName_B <- toupper(net$preferredName_B)
# keep only yellow genes (preferred names)
keep <- net$preferredName_A %in% genes & net$preferredName_B %in% genes
net <- net[keep, ]
write.csv(net, file.path(out_dir, "string_ppi_edges_score400.csv"), row.names = FALSE)
message("PPI edges: ", nrow(net), " unique gene nodes: ",
        length(unique(c(net$preferredName_A, net$preferredName_B))))

# collapse duplicate undirected edges to max score
net$key <- apply(net[, c("preferredName_A", "preferredName_B")], 1, function(x) {
  paste(sort(x), collapse = "|")
})
net <- do.call(rbind, lapply(split(net, net$key), function(d) {
  d[which.max(as.numeric(d$score)), , drop = FALSE]
}))
ends_ab <- strsplit(net$key, "|", fixed = TRUE)
edge_df <- data.frame(
  from = vapply(ends_ab, `[`, character(1), 1),
  to = vapply(ends_ab, `[`, character(1), 2),
  weight = as.numeric(net$score),
  stringsAsFactors = FALSE
)
g <- igraph::graph_from_data_frame(edge_df, directed = FALSE)
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(weight = "max"))
message("Graph: ", igraph::vcount(g), " nodes / ", igraph::ecount(g), " edges")
# ---- MCODE-like vertex weights + seed expansion (Bader & Hogue style) ----
# Parameters aligned with Cytoscape MCODE defaults (approx):
degree_cutoff <- 2
node_score_pct <- 0.2
kcore_min <- 2
fluff <- FALSE

vertex_weight <- function(g, vname) {
  nb <- igraph::neighbors(g, vname)
  if (length(nb) < degree_cutoff) return(0)
  nodes <- unique(c(vname, names(nb)))
  sub <- igraph::induced_subgraph(g, nodes)
  kc <- igraph::coreness(sub)
  kmax <- max(kc)
  if (kmax < kcore_min) return(0)
  core <- igraph::induced_subgraph(sub, igraph::V(sub)$name[kc >= kmax])
  n <- igraph::vcount(core)
  if (n < 2) return(0)
  dens <- igraph::ecount(core) / (n * (n - 1) / 2)
  dens * n
}

vn <- igraph::V(g)$name
vw <- setNames(vapply(vn, function(x) vertex_weight(g, x), numeric(1)), vn)
write.csv(
  data.frame(Gene = names(vw), MCODE_weight = as.numeric(vw), Degree = igraph::degree(g)),
  file.path(out_dir, "mcode_vertex_weights.csv"), row.names = FALSE
)

grow_complex <- function(g, seed, vw, node_score_pct) {
  thr <- vw[seed] * (1 - node_score_pct)
  members <- seed
  frontier <- names(igraph::neighbors(g, seed))
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    cand <- setdiff(frontier, members)
    add <- cand[!is.na(vw[cand]) & vw[cand] >= thr]
    if (length(add)) {
      members <- unique(c(members, add))
      frontier <- unique(c(frontier, unlist(lapply(add, function(x) names(igraph::neighbors(g, x))))))
      changed <- TRUE
    }
  }
  sort(unique(members))
}

# greedily extract non-overlapping complexes from highest-weight seeds
unused <- vn
clusters <- list()
scores <- c()
while (length(unused)) {
  w_sub <- vw[unused]
  if (all(w_sub <= 0) || max(w_sub, na.rm = TRUE) <= 0) break
  seed <- names(which.max(w_sub))
  mem <- grow_complex(g, seed, vw, node_score_pct)
  mem <- intersect(mem, unused)
  if (length(mem) < 3) {
    unused <- setdiff(unused, seed)
    next
  }
  sub <- igraph::induced_subgraph(g, mem)
  # complex score = density * size (MCODE-like report score)
  dens <- if (igraph::vcount(sub) > 1) {
    igraph::ecount(sub) / (igraph::vcount(sub) * (igraph::vcount(sub) - 1) / 2)
  } else {
    0
  }
  sc <- dens * igraph::vcount(sub)
  clusters[[length(clusters) + 1]] <- mem
  scores <- c(scores, sc)
  unused <- setdiff(unused, mem)
}

ord <- order(scores, decreasing = TRUE)
clusters <- clusters[ord]
scores <- scores[ord]
# keep top clusters with score >= 3 (common reporting cutoff) or top 8
keep_n <- max(which(scores >= 3), 0)
if (keep_n < 1) keep_n <- min(5, length(clusters))
keep_n <- min(keep_n, 10, length(clusters))
clusters <- clusters[seq_len(keep_n)]
scores <- scores[seq_len(keep_n)]

clust_df <- do.call(rbind, lapply(seq_along(clusters), function(i) {
  data.frame(
    Cluster = paste0("MCODE", i),
    Score = round(scores[i], 3),
    Size = length(clusters[[i]]),
    Gene = clusters[[i]],
    stringsAsFactors = FALSE
  )
}))
write.csv(clust_df, file.path(out_dir, "mcode_clusters_genes.csv"), row.names = FALSE)

summary_df <- data.frame(
  Cluster = paste0("MCODE", seq_along(clusters)),
  Score = round(scores, 3),
  Size = lengths(clusters),
  Genes = vapply(clusters, function(x) paste(x, collapse = ";"), character(1)),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(out_dir, "mcode_clusters_summary.csv"), row.names = FALSE)
message("MCODE clusters kept: ", nrow(summary_df))
print(summary_df[, c("Cluster", "Score", "Size")])

# ---- Per-cluster enrichment ----
enrich_one <- function(gene_symbols, label) {
  eg <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(eg) < 5) return(list(go = NULL, kegg = NULL))
  go <- tryCatch(
    enrichGO(
      gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
      readable = TRUE
    ),
    error = function(e) NULL
  )
  kegg <- tryCatch(
    enrichKEGG(
      gene = eg$ENTREZID, organism = "hsa",
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2
    ),
    error = function(e) NULL
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  list(go = go, kegg = kegg)
}

# background: all yellow genes mapped
bg_eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

all_go <- list()
all_kegg <- list()
top_terms <- list()
for (i in seq_along(clusters)) {
  lab <- paste0("MCODE", i)
  message("Enrichment: ", lab, " n=", length(clusters[[i]]))
  en <- enrich_one(clusters[[i]], lab)
  if (!is.null(en$go) && nrow(as.data.frame(en$go)) > 0) {
    d <- as.data.frame(en$go)
    d$Cluster <- lab
    all_go[[lab]] <- d
    write.csv(d, file.path(out_dir, paste0(lab, "_GO_BP.csv")), row.names = FALSE)
    top_terms[[lab]] <- head(d$Description, 8)
  } else {
    top_terms[[lab]] <- character(0)
  }
  if (!is.null(en$kegg) && nrow(as.data.frame(en$kegg)) > 0) {
    d <- as.data.frame(en$kegg)
    d$Cluster <- lab
    all_kegg[[lab]] <- d
    write.csv(d, file.path(out_dir, paste0(lab, "_KEGG.csv")), row.names = FALSE)
  }
}

# whole yellow module enrichment for comparison
en_all <- enrich_one(genes, "YELLOW317")
if (!is.null(en_all$go) && nrow(as.data.frame(en_all$go)) > 0) {
  write.csv(as.data.frame(en_all$go), file.path(out_dir, "YELLOW317_GO_BP.csv"), row.names = FALSE)
}
if (!is.null(en_all$kegg) && nrow(as.data.frame(en_all$kegg)) > 0) {
  write.csv(as.data.frame(en_all$kegg), file.path(out_dir, "YELLOW317_KEGG.csv"), row.names = FALSE)
}

# classify clusters: ECM-like vs other
ecm_kw <- c(
  "extracellular matrix", "collagen", "integrin", "cell adhesion",
  "focal adhesion", "actin", "cytoskeleton", "wound healing",
  "TGF", "epithelial to mesenchymal", "ECM"
)
classify_cluster <- function(terms) {
  if (!length(terms)) return("no_sig_enrichment")
  hit <- any(vapply(ecm_kw, function(k) any(grepl(k, terms, ignore.case = TRUE)), logical(1)))
  if (hit) "ECM_adhesion_like" else "possible_nonECM"
}

class_df <- data.frame(
  Cluster = names(top_terms),
  Class = vapply(top_terms, classify_cluster, character(1)),
  TopTerms = vapply(top_terms, function(x) paste(head(x, 5), collapse = " | "), character(1)),
  stringsAsFactors = FALSE
)
class_df <- merge(summary_df[, c("Cluster", "Score", "Size")], class_df, by = "Cluster", all.x = TRUE)
class_df <- class_df[order(-class_df$Score), ]
write.csv(class_df, file.path(out_dir, "mcode_cluster_function_class.csv"), row.names = FALSE)

# ---- Optional ONE-shot prognosis side-check: top 1-3 clusters equal-weight ----
# Reuse step17 style scoring; stop rules same; NOT for iterative fishing.
run_cox_side <- TRUE
if (run_cox_side) {
  suppressPackageStartupMessages(library(survival))
  zmean <- function(mat, gset) {
    gset <- intersect(gset, rownames(mat))
    if (length(gset) < 5) return(rep(NA_real_, ncol(mat)))
    Z <- t(scale(t(mat[gset, , drop = FALSE])))
    colMeans(Z, na.rm = TRUE)
  }
  load_cohort <- function(path, id) {
    x <- readRDS(path)
    list(expr = as.matrix(x$ex), si = x$sampInfo, id = id)
  }
  cohorts <- list(
    load_cohort(file.path(project_root, "data/raw/pdacR/TCGA_PAAD.rds"), "TCGA_PAAD"),
    load_cohort(file.path(project_root, "data/raw/pdacR/Moffitt_GEO_array.rds"), "Moffitt"),
    load_cohort(file.path(project_root, "data/raw/pdacR/Puleo_array.rds"), "Puleo")
  )
  # filter TCGA like before
  tc <- cohorts[[1]]
  si <- tc$si
  keep <- rep(TRUE, nrow(si))
  if ("Filter_by" %in% names(si)) keep <- keep & grepl("whitelist", si$Filter_by, ignore.case = TRUE)
  if ("Histological.type.by.DCC" %in% names(si)) {
    keep <- keep & grepl("ductal", si$Histological.type.by.DCC, ignore.case = TRUE)
  }
  # drop known outlier if present
  oid <- "TCGA-2J-AABV-01A"
  if (oid %in% colnames(tc$expr)) {
    keep[match(oid, colnames(tc$expr))] <- FALSE
  }
  cohorts[[1]]$expr <- tc$expr[, keep, drop = FALSE]
  cohorts[[1]]$si <- si[keep, , drop = FALSE]

  caf_genes <- c("FAP","ACTA2","COL1A1","COL1A2","COL3A1","PDPN","THY1","DCN","LUM","TAGLN","PDGFRB","SPARC")
  n_check <- min(3, length(clusters))
  cox_rows <- list()
  for (i in seq_len(n_check)) {
    gset <- clusters[[i]]
    for (co in cohorts) {
      expr <- co$expr
      rownames(expr) <- toupper(rownames(expr))
      si <- co$si
      sc <- zmean(expr, gset)
      caf <- zmean(expr, caf_genes)
      time <- as.numeric(si$survivalA)
      event <- as.numeric(si$censorA.0yes.1no)
      ok <- is.finite(sc) & is.finite(time) & is.finite(event) & time > 0
      if (sum(ok) < 30) next
      d <- data.frame(time = time[ok], event = event[ok], score = as.numeric(scale(sc[ok])),
                      caf = as.numeric(scale(caf[ok])))
      fit <- coxph(Surv(time, event) ~ score, data = d)
      s <- summary(fit)$coefficients
      ci <- summary(fit)$conf.int
      fit2 <- tryCatch(coxph(Surv(time, event) ~ score + caf, data = d), error = function(e) NULL)
      caf_p <- caf_hr <- NA_real_
      if (!is.null(fit2)) {
        s2 <- summary(fit2)$coefficients
        caf_hr <- exp(s2["score", "coef"])
        caf_p <- s2["score", "Pr(>|z|)"]
      }
      cox_rows[[length(cox_rows) + 1]] <- data.frame(
        Cluster = paste0("MCODE", i),
        Score = scores[i],
        Size = length(gset),
        Cohort = co$id,
        N = nrow(d),
        HR = exp(s[1, "coef"]),
        P = s[1, "Pr(>|z|)"],
        CAF_adj_HR = caf_hr,
        CAF_adj_P = caf_p,
        Rho_vs_CAF = suppressWarnings(cor(d$score, d$caf, use = "complete.obs")),
        stringsAsFactors = FALSE
      )
    }
  }
  cox_all <- do.call(rbind, cox_rows)
  write.csv(cox_all, file.path(out_dir, "mcode_top_clusters_cox_sidecheck.csv"), row.names = FALSE)
}

# ---- Report ----
lines <- c(
  "# Yellow317 MCODE functional dissection",
  paste0("Date: ", as.character(Sys.Date())),
  "",
  "## Framing",
  "- Purpose: internal functional blocks of CAF-linked yellow module (Option 1 supplement).",
  "- NOT a restart of prognostic-signature fishing.",
  "- PPI: STRING physical network, required_score>=400.",
  "- Clustering: MCODE-like vertex weighting + seed expansion (degree_cutoff=2, node_score_pct=0.2, k-core>=2).",
  "",
  "## Clusters (by MCODE-like score)",
  paste(capture.output(print(class_df, row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation rule",
  "- If nearly all high-score clusters are ECM_adhesion_like: confirms CAF-ECM program; use as Supp functional figure.",
  "- If a possible_nonECM cluster appears: discuss as secondary sub-program; still do NOT auto-promote to signature.",
  "- Prognosis side-check: only top 1-3 clusters, equal-weight, one-shot; GO only if external HR>1 P<0.05 AND CAF-adj P<0.05.",
  ""
)
if (exists("cox_all")) {
  lines <- c(
    lines,
    "## Prognosis side-check (top clusters)",
    paste(capture.output(print(cox_all, row.names = FALSE)), collapse = "\n"),
    ""
  )
  # auto verdict
  ext <- cox_all[cox_all$Cohort != "TCGA_PAAD", ]
  ok_go <- any(ext$HR > 1 & ext$P < 0.05 & ext$CAF_adj_P < 0.05, na.rm = TRUE)
  ok_cond <- any(ext$HR > 1 & ext$P < 0.05, na.rm = TRUE) && !ok_go
  verdict <- if (ok_go) {
    "UNEXPECTED GO on a cluster 鈥?freeze list, do not iterate"
  } else if (ok_cond) {
    "CONDITIONAL only 鈥?likely still stroma-linked; keep as functional note"
  } else {
    "NO independent prognosis signal from top MCODE clusters (expected)"
  }
  lines <- c(lines, paste0("## Side-check verdict: ", verdict), "")
}
lines <- c(
  lines,
  "## Bottom line",
  "Use MCODE blocks to answer: what functional subdivisions exist inside yellow CAF-ECM?",
  "Do not use MCODE as a gene-picker for a new OS model unless a non-ECM cluster clearly passes the locked side-check."
)
writeLines(lines, file.path(out_dir, "STEP18_YELLOW317_MCODE_CN.md"))
message("Done. Output: ", out_dir)
