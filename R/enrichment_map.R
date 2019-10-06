#' Title
#'
#' TODO
#'
#' @param res_enrich TODO
#' @param res_de TODO
#' @param annotation_obj TODO
#' @param n_gs TODO
#' @param overlap_threshold TODO
#' @param scale_edges_width TODO
#' @param color_by TODO
#' @param size_by TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#' @param genesetid_colname TODO
#'
#' TODOTODO: similarity measures, say, jaccard, or simple overlap
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
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
  enriched_gsdescs <- vapply(enriched_gsids, function(arg) Definition(GOTERM[[arg]]))

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
        overlap_jaccard_index(unlist(enrich2list[enriched_gsids[i]]),
                              unlist(enrich2list[enriched_gsids[j]]))
    }
  }

  # oooooor TODOTODO: go directly from this adjacency matrix? use overlap as weights
  om_df <- as.data.frame(overlap_matrix)
  om_df$id <- rownames(om_df)

  omm <- pivot_longer(om_df, seq_len(n_gs))
  colnames(omm) <- c("gs_1","gs_2","value")
  # eliminate rows of diagonal...
  omm <- omm[omm$gs_1 != omm$gs_2, ]
  # ... and the ones from the other triangular portion
  omm <- omm[!is.na(omm$value), ]

  # omm <- reshape2::melt(overlap_matrix)
  # omm <- omm[omm$Var1 != omm$Var2, ]
  # omm <- omm[!is.na(omm$value), ]

  # use this to construct the graph
  g <- graph.data.frame(omm[,c(1,2)], directed = FALSE)

  E(g)$width <- sqrt(omm$value * scale_edges_width)
  g <- delete.edges(g, E(g)[omm$value < overlap_threshold])

  idx <- match(V(g)$name, res_enrich$Term)

  gs_size <- vapply(enrich2list[idx], length)

  V(g)$size <- 5 * sqrt(gs_size)
  V(g)$original_size <- gs_size


  colVar <- res_enrich[idx, color_by]
  if(all(colVar < 1)) # likely p-values...
    colVar <- -log10(colVar)
  # V(g)$color <- colVar

  mypal <- (scales::alpha(
    colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu",11))(50), 0.8))
  mypal_hover <- (scales::alpha(
    colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu",11))(50), 0.5))
  mypal_select <- (scales::alpha(
    colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu",11))(50), 1))

  # V(g)$color <- map2color(colVar,mypal,limits = range(colVar))
  V(g)$color.background <- map2color(colVar,mypal,limits = range(colVar))
  V(g)$color.highlight <- map2color(colVar,mypal_select,limits = range(colVar))
  V(g)$color.hover <- map2color(colVar,mypal_hover,limits = range(colVar))

  # TODOTODO: some kind of border prettifying?
  V(g)$color.border <- "black"

  # additional specification of edge colors
  E(g)$color <- "lightgrey"

  # TODOTODO: think of sorting nodes alphabetically? ..
  # g <- permute.vertices(g,Matrix::invPerm(order(V(g)$name)))
  return(g)
}

# g %>% visIgraph() %>% visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = TRUE)



