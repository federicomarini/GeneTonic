enrichment_map <- function(res_enrich,
                           res_de,
                           annotation_obj,
                           n_gs = 50,
                           overlap_threshold = 0.1,
                           scale_edges_width = 200,
                           color_by = "p.value_elim",
                           size_by,
                           genes_colname = "genes",
                           genesetname_colname = "Term",
                           genesetid_colname = "GO.ID") {

  # if we want to allow for other feats to be colored by, check that some aggregated scores are there
  # TODOTODO, if... otherwise, compute aggr scores

  enriched_gsids <- res_enrich[[genesetid_colname]]
  enriched_gsnames <- res_enrich[[genesetname_colname]]
  enriched_gsdescs <- sapply(enriched_gsids, function(arg) Definition(GOTERM[[arg]]))

  rownames(res_enrich) <- enriched_gsids

  enrich2list <- lapply(seq_len(n_gs),function(gs){
    # goterm <- res_enrich$Term[gs]
    go_genes <- res_enrich$genes[gs]
    go_genes <- strsplit(go_genes,",") %>% unlist
    return(go_genes)
  })
  names(enrich2list) <- enriched_gsids[seq_len(n_gs)]

  n <- length(enrich2list)
  overlap_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(overlap_matrix) <- colnames(overlap_matrix) <- enriched_gsnames[seq_len(n_gs)]

  for (i in 1:n) {
    # no need to work on full mat, it is simmetric
    for (j in i:n) {
      overlap_matrix[i, j] <-
        overlap_ratio(enrich2list[enriched_gsids[i]], enrich2list[enriched_gsids[j]])
    }
  }

  omm <- reshape2::melt(overlap_matrix)
  # eliminate rows of diagonal...
  omm <- omm[omm$Var1 != omm$Var2, ]
  # ... and the ones from the other triangular portion
  omm <- omm[!is.na(omm$value), ]

  # use this to construct the graph
  g <- graph.data.frame(omm[,c(1,2)], directed = FALSE)

  E(g)$width <- sqrt(omm$value * scale_edges_width)
  g <- delete.edges(g, E(g)[omm$value < overlap_threshold])

  idx <- match(V(g)$name, res_enrich$Term)

  gs_size <- sapply(enrich2list[idx], length)

  V(g)$size <- 5 * sqrt(gs_size)
  V(g)$original_size <- gs_size


  colVar <- res_enrich[idx, color_by]
  if(all(colVar < 1)) # likely p-values...
    colVar <- -log10(colVar)
  V(g)$color <- colVar

  mypal <- (scales::alpha(
    colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu",11))(50), 0.8))
  V(g)$color <- map2color(colVar,mypal,limits = range(colVar))


  # additional specification of edge colors
  E(g)$color <- "lightgrey"

  return(g)
}

# g %>% visIgraph() %>%
  # visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = TRUE)

overlap_ratio <- function (x, y)
{
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x, y)))
}


map2color <- function(x, pal, limits = NULL) {
  if(is.null(limits))
    limits=range(x)
  pal[findInterval(x, seq(limits[1],
                          limits[2],
                          length.out=length(pal)+1),
                   all.inside=TRUE)]
}


