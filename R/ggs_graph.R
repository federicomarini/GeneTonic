#' Construct a gene-geneset-graph
#'
#' Construct a gene-geneset-graph from the results of a functional enrichment
#' analysis
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to see the
#' formatting requirements.
#' @param res_de A `DESeqResults` object.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be included
#' @param genes_colname Character, specifying which column of the `res_enrich`
#' object contains the genes assigned to each gene set, detected as differentially
#' expressed. Defaults to `genes`.
#' @param genesetname_colname Character, specifies which column of the `res_enrich`
#' object contains a description of the gene set. Defaults to `Term`.
#' @param genesetid_colname Character, specifies which column of the `res_enrich`
#' object contains a unique identifier of the gene set. Defaults to `GO.ID`.
#' @param prettify Logical, controlling the aspect of the returned graph object.
#' If TRUE (default value), different shapes of the nodes are returned, based on
#' the node type
#' @param geneset_graph_color Character value, specifying which color should be
#' used for the fill of the shapes related to the gene sets.
#' @param genes_graph_colpal TODO
#'
#' @return An `igraph` obejct to be further manipulated or processed/plotted (e.g.
#' via [igraph::plot.igraph()] or [visNetwork::visIgraph()])
#' @export
#'
#' @examples
#' #TODO
# TODO: maybe rename to ggs_network? I like it so!
ggs_graph <- function(res_enrich,
                      res_de,
                      annotation_obj = NULL,
                      n_gs = 15,
                      genes_colname = "genes",
                      genesetname_colname = "Term",
                      genesetid_colname = "GO.ID",
                      prettify = TRUE,
                      geneset_graph_color = "gold",
                      genes_graph_colpal) {

  # res_enrich has to have a column called containing the genes annotated to the term
  # TODOTODO

  # verify the genesets are sorted in a meaningful way?
  #TODOTODO

  enriched_gsids <- res_enrich[[genesetid_colname]]
  enriched_gsnames <- res_enrich[[genesetname_colname]]
  enriched_gsdescs <- vapply(enriched_gsids,
                             function(arg) Definition(GOTERM[[arg]]),
                             character(1))

  enrich2list <- lapply(seq_len(n_gs), function(gs) {
    # goterm <- res_enrich$Term[gs]
    go_genes <- res_enrich$genes[gs]
    go_genes <- strsplit(go_genes, ",") %>% unlist
    return(go_genes)
  })
  names(enrich2list) <- enriched_gsnames[seq_len(n_gs)]

  list2df <- lapply(seq_len(length(enrich2list)), function(gs) {
    data.frame(
      gsid = rep(names(enrich2list[gs]), length(enrich2list[[gs]])),
      gene = enrich2list[[gs]])
  })
  list2df <- do.call("rbind", list2df)

  g <- graph.data.frame(list2df, directed = FALSE)

  nodeIDs_gs <- which(names(V(g)) %in% enriched_gsnames)
  nodeIDs_genes <- which(!(names(V(g)) %in% enriched_gsnames))


  if (prettify) {
    # different shapes based on the node type

    V(g)$nodetype <- NA
    V(g)$nodetype[nodeIDs_gs] <- "GeneSet"
    V(g)$nodetype[nodeIDs_genes] <- "Feature"

    # TODOTODO: this does not work with visNetwork?
    # V(g)$value <- 15 # size? size2? or does this not work with the shapes I selected?
    # V(g)$value[nodeIDs_gs] <- 45

    # V(g)$goid <- NA
    # V(g)$goid[1:n] <- goids

    # this one is handled correctly by visNetwork
    V(g)$shape <- c("box", "ellipse")[factor(V(g)$nodetype, levels = c("GeneSet", "Feature"))]


    # different colors for the gene nodes f logFC
    fcs_genes <- res_de[annotation_obj$gene_id[match((V(g)$name[nodeIDs_genes]), annotation_obj$gene_name)], ]$log2FoldChange


    mypal <- rev(scales::alpha(
      colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.4))
    V(g)$color[nodeIDs_genes] <- map2color(fcs_genes, mypal, limits = c(-4, 4))
    V(g)$color[nodeIDs_gs] <- geneset_graph_color

    # title for tooltips
    V(g)$title <- NA
    V(g)$title[nodeIDs_gs] <- paste0("<h4>",
                                     sprintf('<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank">%s</a>', enriched_gsids[nodeIDs_gs], enriched_gsids[nodeIDs_gs]), "</h4><br>",
                                     V(g)$name[nodeIDs_gs], "<br><br>",
                                     enriched_gsdescs[nodeIDs_gs])
    V(g)$title[nodeIDs_genes] <- paste0("<h4>", V(g)$name[nodeIDs_genes], "</h4><br>",
                                        "logFC = ", format(round(fcs_genes, 2), nsmall = 2))


  } else {

    V(g)$color[nodeIDs_genes] <- "#B3B3B3"
    V(g)$color[nodeIDs_gs] <- "#E5C494"
  }
  return(g)
}
