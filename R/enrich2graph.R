#' Title
#'
#' TODO
#'
#' @param res_enrich TODO
#' @param res_de TODO
#' @param annotation_obj TODO
#' @param n_gs TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#' @param genesetid_colname TODO
#' @param prettify TODO
#' @param geneset_graph_color TODO
#' @param genes_graph_colpal TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' #TODO
# TODO: maybe rename to ggs_network? I like it so!
enrich2graph <- function(res_enrich,
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
