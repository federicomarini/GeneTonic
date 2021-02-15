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
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be included
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be included in addition to 
#' the top ones (via `n_gs`)
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
                      gtl = NULL,
                      n_gs = 15,
                      gs_ids = NULL,
                      prettify = TRUE,
                      geneset_graph_color = "gold",
                      genes_graph_colpal = NULL) {
  
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }
  
  stopifnot(is.numeric(n_gs))
  stopifnot(is.logical(prettify))
  
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
                             function(arg) tryCatch(
                               Definition(GOTERM[[arg]]),
                               error = function(e) "--- GO Term not found ---"
                               ),
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







#' Extract the backbone for the gene-geneset graph
#' 
#' Extract the backbone for the gene-geneset graph, either for the genes or for the 
#' genesets
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param res_de A `DESeqResults` object.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be included
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be included in addition to 
#' the top ones (via `n_gs`)
#' @param bb_on A character string, either "genesets" or "features", to specify which
#' entity should be based the backbone graph on.
#' @param bb_method A character string, referring to the function to be called (
#' from the `backbone` package) for computing the backbone of the specified 
#' bipartite graph. Defaults to "sdsm", as recommended in the `backbone` package.
#' @param bb_extract_alpha A numeric value, specifying the significance level to 
#' use when detecting the backbone of the network 
#' @param bb_extract_fwer A character string, defaulting to "none", specifying 
#' which method to use for the multiple testing correction for controlling the 
#' family-wise error rate
#' @param bb_fullinfo Logical value, determining what will be returned as output:
#' either a simple `ìgraph` object with the graph backbone (if set to `FALSE`), 
#' or a list object containing also the `backbone` object, and the gene-geneset 
#' graph used for the computation (if `TRUE`)
#' @param bb_remove_singletons Logical value, defines whether to remove or leave 
#' in the returned graph the nodes that are not connected to other vertices
#' @param color_graph Logical value, specifies whether to use information about 
#' genesets or features to colorize the nodes, e.g. for this info to be used in 
#' interactive versions of the graph
#' @param color_by_geneset Character string, corresponding to the column in 
#' `res_enrich` to be used for coloring the nodes if `bb_on` is set to "genesets".
#' Defaults to the "z_score", which can be obtained via `get_aggrscores()`
#' @param color_by_feature Character string, corresponding to the column in 
#' `res_de` to be used for coloring the nodes if `bb_on` is set to "features".
#' Defaults to the "log2FoldChange", which should be normally included in a
#' DESeqResults object.
#' @param ... Additional parameters to be passed internally
#'
#' @return According to the `bb_fullinfo`, either a simple `ìgraph` object with 
#' the graph backbone, or a named list object containing:
#' - the `igraph` of the extracted backbone
#' - the `backbone` object itself
#' - the gene-geneset graph used for the computation
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
#' ggs_bbg <- ggs_backbone(res_enrich,
#'                              res_de,
#'                              anno_df,
#'                              n_gs = 50,
#'                              bb_on = "genesets",
#'                              color_graph = TRUE,
#'                              color_by_geneset = "z_score"
#' )
#' plot(ggs_bbg)
#' 
#' # if desired, one can also plot the interactive version
#' visNetwork::visIgraph(ggs_bbg)
#' 
ggs_backbone <- function(res_enrich,
                         res_de,
                         annotation_obj = NULL,
                         gtl = NULL,
                         n_gs = 15,
                         gs_ids = NULL,
                         bb_on = c("genesets", "features"),
                         bb_method = c("sdsm", "fsdm", "hyperg"),
                         bb_extract_alpha = 0.05,
                         bb_extract_fwer = c("none","bonferroni","holm"),
                         bb_fullinfo = FALSE,
                         bb_remove_singletons = TRUE,
                         color_graph = TRUE,
                         color_by_geneset = "z_score",
                         color_by_feature = "log2FoldChange",
                         ...) {
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }
  
  stopifnot(is.logical(bb_fullinfo))
  stopifnot(is.logical(bb_remove_singletons))
  stopifnot(is.logical(color_graph))
  
  # check that columns to encode the colors are present
  if (color_graph) {
    if (bb_on == "genesets") {
      if (!color_by_geneset %in% colnames(res_enrich))
        stop("Your res_enrich object does not contain the ",
             color_by_geneset,
             " column.\n",
             "Compute this first or select another column to use for the color.")
      
    } else if (bb_on == "features") {
      if (!color_by_feature %in% colnames(res_de))
        stop("Your res_de object does not contain the ",
             color_by_feature,
             " column.\n",
             "Compute this first or select another column to use for the color.")
    }
  }
  
  bb_method <- match.arg(bb_method, c("sdsm", "fsdm", "hyperg"))
  bb_extract_fwer = match.arg(bb_extract_fwer, c("none","bonferroni","holm"))
  bb_on <- match.arg(bb_on, c("genesets", "features"))
  
  # first, compute the ggs graph object
  ggs <- ggs_graph(res_enrich = res_enrich,
                   res_de = res_de,
                   annotation_obj = annotation_obj,
                   gtl = gtl,
                   n_gs = n_gs,
                   gs_ids = gs_ids)
  
  # for making this a formal bipartite graph
  V(ggs)$type <- V(ggs)$nodetype=="GeneSet"
  
  bpm <- igraph::as_incidence_matrix(ggs)
  
  if (bb_on =="features") {
    bpm_for_backbone <- bpm
  } else if (bb_on == "genesets") {
    bpm_for_backbone <- t(bpm)
  }
  
  if (bb_method == "sdsm") {
    bbobj <- backbone::sdsm(bpm_for_backbone)
  } else if (bb_method == "fdsm") {
    bbobj <- backbone::fdsm(bpm_for_backbone, trials = 1000)
  } else if (bb_method == "hyperg") {
    bbobj <- backbone::hyperg(bpm_for_backbone)
  }
  
  bbextracted <- backbone::backbone.extract(bbobj,
                                            alpha = bb_extract_alpha,
                                            fwer = bb_extract_fwer)
  
  bbgraph <- igraph::graph_from_adjacency_matrix(bbextracted, mode = "undirected")
  
  if (bb_remove_singletons) {
    bbgraph <- igraph::delete_vertices(bbgraph, !(igraph::degree(bbgraph) >= 1 ))
  }
  
  if (color_graph) {
    if (bb_on == "genesets") {
      # will use the summarized Z-score
      color_by <- color_by_geneset
      idx <- match(V(bbgraph)$name, res_enrich$gs_description)
      col_var <- res_enrich[idx, color_by]
      
      mypal <- rev(scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.8))
      mypal_hover <- rev(scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.5))
      mypal_select <- rev(scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 1))
      
      V(bbgraph)$color.background <- map2color(col_var, mypal, limits = range(col_var))
      V(bbgraph)$color.highlight <- map2color(col_var, mypal_select, limits = range(col_var))
      V(bbgraph)$color.hover <- map2color(col_var, mypal_hover, limits = range(col_var))
      
      V(bbgraph)$color.border <- "black"
      
      # additional specification of edge colors
      E(bbgraph)$color <- "lightgrey"
      
    } else if (bb_on == "features") {
      # will use the logFoldChange
      color_by <- color_by_feature
      idx <- annotation_obj$gene_id[match(V(bbgraph)$name, annotation_obj$gene_name)]
      col_var <- res_de[idx, color_by]
      
      mypal <- rev(scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.8))
      mypal_hover <- rev(scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.5))
      mypal_select <- rev(scales::alpha(
        colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 1))
      
      V(bbgraph)$color.background <- map2color(col_var, mypal, limits = range(col_var))
      V(bbgraph)$color.highlight <- map2color(col_var, mypal_select, limits = range(col_var))
      V(bbgraph)$color.hover <- map2color(col_var, mypal_hover, limits = range(col_var))
      
      V(bbgraph)$color.border <- "black"
      
      # additional specification of edge colors
      E(bbgraph)$color <- "lightgrey"
    }
    
  }
  
  if (bb_fullinfo) {
    return(
      list(
        bbgraph = bbgraph,
        bbobj = bbobj,
        ggs = ggs
      )
    )
  } else {
    return(bbgraph)
  }
}  
