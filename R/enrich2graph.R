#' Title
#'
#' @param res_enrich TODO
#' @param res_de TODO
#' @param n_nodes TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#' @param genesetid_colname TODO
#' @param prettify TODO
#' @param geneset_graph_color TODO
#' @param genes_graph_colpal TODO
#' @param annotation_obj TODO
#'
#' @return
#' @export
#'
#' @examples
enrich2graph <- function(res_enrich,
                         res_de,
                         n_nodes = 15,
                         genes_colname = "genes",
                         genesetname_colname = "Term",
                         genesetid_colname = "GO.ID",
                         prettify = TRUE,
                         geneset_graph_color = "gold",
                         genes_graph_colpal,
                         annotation_obj = NULL
                         ) {

  # res_enrich has to have a column called containing the genes annotated to the term
  # TODOTODO

  # verify the genesets are sorted in a meaningful way?
  #TODOTODO

  enriched_gsids <- res_enrich[[genesetid_colname]]
  enriched_gsnames <- res_enrich[[genesetname_colname]]

  enrich2list <- lapply(seq_len(n_nodes),function(gs){
    # goterm <- res_enrich$Term[gs]
    go_genes <- res_enrich$genes[gs]
    go_genes <- strsplit(go_genes,",") %>% unlist
    return(go_genes)
  })
  names(enrich2list) <- enriched_gsnames[seq_len(n_nodes)]

  list2df <- lapply(seq_len(length(enrich2list)), function(gs) {
    data.frame(
      gsid = rep(names(enrich2list[gs]), length(enrich2list[[gs]])),
      gene = enrich2list[[gs]])
  })
  list2df <- do.call("rbind", list2df)

  g <- graph.data.frame(list2df, directed = FALSE)

  if(prettify){
    # different shapes based on the node type
    nodeIDs_gs <- which(names(V(g)) %in% enriched_gsnames)
    nodeIDs_genes <- which(!(names(V(g)) %in% enriched_gsnames))

    V(g)$nodetype <- NA
    V(g)$nodetype[nodeIDs_gs] <- "GeneSet"
    V(g)$nodetype[nodeIDs_genes] <- "Feature"

    # V(g)$goid <- NA
    # V(g)$goid[1:n] <- goids

    # this one is handled correctly by visNetwork
    V(g)$shape <- c("box","ellipse")[factor(V(g)$nodetype, levels=c("GeneSet","Feature"))]


    # different colors for the gene nodes f logFC
    fcs_genes <- res_de$log2FoldChange[match((V(g)$name[nodeIDs_genes]),annotation_obj$gene_name)]

    map2color <- function(x, pal, limits = NULL) {
      if(is.null(limits))
        limits=range(x)
      pal[findInterval(x, seq(limits[1],
                              limits[2],
                              length.out=length(pal)+1),
                       all.inside=TRUE)]
    }
    mypal <- rev(scales::alpha(
      colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu",11))(50), 0.4))
    V(g)$color[nodeIDs_genes] <- map2color(fcs_genes,mypal,limits = c(-4,4))
    V(g)$color[nodeIDs_gs] <- geneset_graph_color
    } else {

      V(g)$color[nodeIDs_genes] <- "#B3B3B3"
      V(g)$color[nodeIDs_gs] <- "#E5C494"
  }
  return(g)

}

