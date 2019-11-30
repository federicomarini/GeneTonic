#' Creates an enrichment map for the results of functional enrichment
#'
#' Generates a graph for the enrichment map, combining information from `res_enrich`
#' and `res_de`. This object can be further plotted, e.g. statically via
#' [igraph::plot.igraph()], or dynamically via [visNetwork::visIgraph()]
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param res_de A `DESeqResults` object.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be displayed
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be displayed.
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
#' library("macrophage")
#' library("DESeq2")
#' library("org.Hs.eg.db")
#' library("AnnotationDbi")
#'
#' # dds object
#' data("gse", package = "macrophage")
#' dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#' # annotation object
#' anno_df <- data.frame(
#'   gene_id = rownames(dds_macrophage),
#'   gene_name = mapIds(org.Hs.eg.db,
#'                      keys = rownames(dds_macrophage),
#'                      column = "SYMBOL",
#'                      keytype = "ENSEMBL"),
#'   stringsAsFactors = FALSE,
#'   row.names = rownames(dds_macrophage)
#' )
#'
#' # res object
#' data(res_de_macrophage, package = "GeneTonic")
#' res_de <- res_macrophage_IFNg_vs_naive
#'
#' # res_enrich object
#' data(res_enrich_macrophage, package = "GeneTonic")
#' res_enrich <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
#' res_enrich <- get_aggrscores(res_enrich, res_de, anno_df)
#'
#' em <- enrichment_map(res_enrich,
#'                      res_de,
#'                      anno_df,
#'                      n_gs = 20
#' )
#'
#' em
#'
#' # could be viewed interactively with
#' # library(visNetwork)
#' # library(magrittr)
#' # em %>%
#' #   visIgraph() %>%
#' #   visOptions(highlightNearest = list(enabled = TRUE,
#' #                                      degree = 1,
#' #                                      hover = TRUE),
#' #             nodesIdSelection = TRUE)
enrichment_map <- function(res_enrich,
                           res_de,
                           annotation_obj,
                           n_gs = 50,
                           gs_ids = NULL,
                           overlap_threshold = 0.1,
                           scale_edges_width = 200,
                           color_by = "gs_pvalue",
                           size_by) {

  # if we want to allow for other feats to be colored by, check that some aggregated scores are there
  # TODOTODO, if... otherwise, compute aggr scores

  # enriched_gsids <- res_enrich$gs_id
  # enriched_gsnames <- res_enrich$gs_description
  # enriched_gsdescs <- vapply(enriched_gsids, function(arg) Definition(GOTERM[[arg]]), character(1))

  n_gs <- min(n_gs, nrow(res_enrich))

  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)],  # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id]  # the ones specified from the custom list
    )
  )

  overlap_matrix <- create_jaccard_matrix(res_enrich,
                                          n_gs = n_gs,
                                          gs_ids = gs_ids,
                                          return_sym = FALSE)

  rownames(overlap_matrix) <- colnames(overlap_matrix) <- res_enrich[rownames(overlap_matrix), "gs_description"]

  # oooooor TODOTODO: go directly from this adjacency matrix? use overlap as weights
  om_df <- as.data.frame(overlap_matrix)
  om_df$id <- rownames(om_df)

  omm <- pivot_longer(om_df, seq_len(length(gs_to_use)))
  colnames(omm) <- c("gs_1", "gs_2", "value")
  # eliminate rows of diagonal...
  omm <- omm[omm$gs_1 != omm$gs_2, ]
  # ... and the ones from the other triangular portion
  omm <- omm[!is.na(omm$value), ]

  # omm <- reshape2::melt(overlap_matrix)
  # omm <- omm[omm$Var1 != omm$Var2, ]
  # omm <- omm[!is.na(omm$value), ]

  # use this to construct the graph
  emg <- graph.data.frame(omm[, c(1,2)], directed = FALSE)

  E(emg)$width <- sqrt(omm$value * scale_edges_width)
  emg <- delete.edges(emg, E(emg)[omm$value < overlap_threshold])

  idx <- match(V(emg)$name, res_enrich$gs_description)

  gs_size <- res_enrich$gs_de_count[idx]

  V(emg)$size <- 5 * sqrt(gs_size)
  V(emg)$original_size <- gs_size


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
  V(emg)$color.background <- map2color(col_var, mypal, limits = range(col_var))
  V(emg)$color.highlight <- map2color(col_var, mypal_select, limits = range(col_var))
  V(emg)$color.hover <- map2color(col_var, mypal_hover, limits = range(col_var))

  # TODOTODO: some kind of border prettifying?
  V(emg)$color.border <- "black"

  # additional specification of edge colors
  E(emg)$color <- "lightgrey"

  # re-sorting the vertices alphabetically
  rank_gs <- rank(V(emg)$name)
  emg <- permute.vertices(emg, rank_gs)

  return(emg)
}

# g %>% visIgraph() %>% visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = TRUE)
