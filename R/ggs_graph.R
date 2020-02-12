#' Construct a gene-geneset-graph
#'
#' Construct a gene-geneset-graph from the results of a functional enrichment
#' analysis
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param res_de A `DESeqResults` object.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be included
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be displayed.
#' @param prettify Logical, controlling the aspect of the returned graph object.
#' If TRUE (default value), different shapes of the nodes are returned, based on
#' the node type
#' @param geneset_graph_color Character value, specifying which color should be
#' used for the fill of the shapes related to the gene sets.
#' @param genes_graph_colpal A vector of colors, also provided with their hex
#' string, to be used as a palette for coloring the gene nodes. If unspecified,
#' defaults to a color ramp palette interpolating from blue through yellow to red.
#'
#' @return An `igraph` object to be further manipulated or processed/plotted (e.g.
#' via [igraph::plot.igraph()] or 
#' [visNetwork::visIgraph()][visNetwork::visNetwork-igraph])
#' @export
#'
#' @examples
#'
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
#' ggs <- ggs_graph(res_enrich,
#'                  res_de,
#'                  anno_df
#'                 )
#'
#' ggs
#'
#' #' # could be viewed interactively with
#' # library(visNetwork)
#' # library(magrittr)
#' # ggs %>%
#' #   visIgraph() %>%
#' #   visOptions(highlightNearest = list(enabled = TRUE,
#' #                                      degree = 1,
#' #                                      hover = TRUE),
#' #             nodesIdSelection = TRUE)
ggs_graph <- function(res_enrich,
                      res_de,
                      annotation_obj = NULL,
                      n_gs = 15,
                      gs_ids = NULL,
                      prettify = TRUE,
                      geneset_graph_color = "gold",
                      genes_graph_colpal = NULL) {

  if (!is.null(genes_graph_colpal)) {

    if (!is(genes_graph_colpal, "character"))
      stop("Please check that you are correctly providing the color palette, ",
           "it should be encoded as a vector of colors specified as characters ",
           "(textual or hex codes)")

    if (!all(check_colors(genes_graph_colpal)))
      stop("You are providing your color palette in a format which ",
           "\ncan not be handled by `grDevices::col2rgb`. \n\n",
           "Try running `check_colors` on the palette object.")
  }

  n_gs <- min(n_gs, nrow(res_enrich))

  enriched_gsids <- res_enrich[["gs_id"]]
  enriched_gsnames <- res_enrich[["gs_description"]]
  enriched_gsdescs <- vapply(enriched_gsids,
                             function(arg) Definition(GOTERM[[arg]]),
                             character(1))

  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)],  # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id]  # the ones specified from the custom list
    )
  )

  enrich2list <- lapply(gs_to_use, function(gs) {
    go_genes <- res_enrich[gs, "gs_genes"]
    go_genes <- unlist(strsplit(go_genes, ","))
    return(go_genes)
  })
  names(enrich2list) <- res_enrich[gs_to_use, "gs_description"]

  list2df <- lapply(seq_along(enrich2list), function(gs) {
    data.frame(
      gsid = rep(names(enrich2list[gs]), length(enrich2list[[gs]])),
      gene = enrich2list[[gs]])
  })
  list2df <- do.call("rbind", list2df)

  g <- graph.data.frame(list2df, directed = FALSE)

  nodeIDs_gs <- which(names(V(g)) %in% enriched_gsnames)
  nodeIDs_genes <- which(!(names(V(g)) %in% enriched_gsnames))

  V(g)$nodetype <- NA
  V(g)$nodetype[nodeIDs_gs] <- "GeneSet"
  V(g)$nodetype[nodeIDs_genes] <- "Feature"

  if (prettify) {
    # different shapes based on the node type

    # this does not work with visNetwork?
    # V(g)$value <- 15 # size? size2? or does this not work with the shapes I selected?
    # V(g)$value[nodeIDs_gs] <- 45

    # this one is handled correctly by visNetwork
    V(g)$shape <- c("box", "ellipse")[factor(V(g)$nodetype, levels = c("GeneSet", "Feature"))]


    # different colors for the gene nodes in function of their logFC
    fcs_genes <- res_de[annotation_obj$gene_id[match((V(g)$name[nodeIDs_genes]), annotation_obj$gene_name)], ]$log2FoldChange

    if (!is.null(genes_graph_colpal)) {
      mypal <- genes_graph_colpal
    } else {
      mypal <- rev(scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.4))
    }

    V(g)$color[nodeIDs_genes] <- map2color(fcs_genes, mypal, limits = c(-4, 4))
    V(g)$color[nodeIDs_gs] <- geneset_graph_color

    # title for tooltips
    V(g)$title <- NA
    V(g)$title[nodeIDs_gs] <- paste0("<h4>",
                                     sprintf('<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank">%s</a>', enriched_gsids[nodeIDs_gs], enriched_gsids[nodeIDs_gs]), "</h4><br>",
                                     V(g)$name[nodeIDs_gs], "<br><br>",
                                     paste0(strwrap(enriched_gsdescs[nodeIDs_gs], 50), collapse = "<br>"))
    V(g)$title[nodeIDs_genes] <- paste0("<h4>", V(g)$name[nodeIDs_genes], "</h4><br>",
                                        "logFC = ", format(round(fcs_genes, 2), nsmall = 2))


  } else {

    V(g)$color[nodeIDs_genes] <- "#B3B3B3"
    V(g)$color[nodeIDs_gs] <- "#E5C494"
  }

  # re-sorting the vertices alphabetically
  rank_gs <- rank(V(g)$name[V(g)$nodetype == "GeneSet"])
  rank_feats <- rank(V(g)$name[V(g)$nodetype == "Feature"]) +
    length(rank_gs) # to keep the GeneSets first
  g <- permute.vertices(g, c(rank_gs, rank_feats))

  return(g)
}
