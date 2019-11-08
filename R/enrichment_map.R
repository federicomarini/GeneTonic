#' Creates an enrichment map for the results of functional enrichment
#'
#' Generates a graph for the enrichment map, combining information from `res_enrich`
#' and `res_de`. This object can be further plotted, e.g. statically via
#' [igraph::plot.igraph()], or dynamically via [visNetwork::visIgraph()]
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to see the
#' formatting requirements
#' @param res_de A `DESeqResults` object.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be displayed
#' @param overlap_threshold Numeric value, between 0 and 1. Defines the threshold
#' to be used for removing edges in the enrichment map - edges below this value
#' will be excluded from the final graph. Defaults to 0.1.
#' @param scale_edges_width TODO
#' @param color_by TODO
#' @param size_by TODO
#'
#' TODOTODO: similarity measures, say, jaccard, or simple overlap
#'
#' @return An `igraph` object to be further manipulated or processed/plotted
#'
#' @seealso [GeneTonic()] embeds an interactive visualization for the enrichment map
#'
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
                           color_by = "gs_pvalue",
                           size_by) {

  # if we want to allow for other feats to be colored by, check that some aggregated scores are there
  # TODOTODO, if... otherwise, compute aggr scores

  enriched_gsids <- res_enrich$gs_id
  enriched_gsnames <- res_enrich$gs_description
  enriched_gsdescs <- vapply(enriched_gsids, function(arg) Definition(GOTERM[[arg]]), character(1))

  rownames(res_enrich) <- enriched_gsids

  enrich2list <- lapply(seq_len(n_gs), function(gs) {
    go_genes <- res_enrich$gs_genes[gs]
    go_genes <- strsplit(go_genes, ",") %>% unlist
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
  colnames(omm) <- c("gs_1", "gs_2", "value")
  # eliminate rows of diagonal...
  omm <- omm[omm$gs_1 != omm$gs_2, ]
  # ... and the ones from the other triangular portion
  omm <- omm[!is.na(omm$value), ]

  # omm <- reshape2::melt(overlap_matrix)
  # omm <- omm[omm$Var1 != omm$Var2, ]
  # omm <- omm[!is.na(omm$value), ]

  # use this to construct the graph
  g <- graph.data.frame(omm[, c(1,2)], directed = FALSE)

  E(g)$width <- sqrt(omm$value * scale_edges_width)
  g <- delete.edges(g, E(g)[omm$value < overlap_threshold])

  idx <- match(V(g)$name, res_enrich$gs_description)

  gs_size <- vapply(enrich2list[idx], length, numeric(1))

  V(g)$size <- 5 * sqrt(gs_size)
  V(g)$original_size <- gs_size


  col_var <- res_enrich[idx, color_by]
  if (all(col_var < 1)) # likely p-values...
    col_var <- -log10(col_var)
  # V(g)$color <- colVar

  # TODO: palette changes if it is z_score VS pvalue
  mypal <- (scales::alpha(
    colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.8))
  mypal_hover <- (scales::alpha(
    colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.5))
  mypal_select <- (scales::alpha(
    colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 1))

  # V(g)$color <- map2color(colVar,mypal,limits = range(colVar))
  V(g)$color.background <- map2color(col_var, mypal, limits = range(col_var))
  V(g)$color.highlight <- map2color(col_var, mypal_select, limits = range(col_var))
  V(g)$color.hover <- map2color(col_var, mypal_hover, limits = range(col_var))

  # TODOTODO: some kind of border prettifying?
  V(g)$color.border <- "black"

  # additional specification of edge colors
  E(g)$color <- "lightgrey"

  # TODOTODO: think of sorting nodes alphabetically? ..
  # g <- permute.vertices(g,Matrix::invPerm(order(V(g)$name)))
  ## TODO: or simply rank() :)
  return(g)
}

# g %>% visIgraph() %>% visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = TRUE)
