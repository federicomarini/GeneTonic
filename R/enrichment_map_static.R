#' Creates an enrichment map for the results of functional enrichment
#'
#' Generates a graph for the enrichment map, combining information from `res_enrich`
#' and `res_de`. This object can be further plotted, e.g. statically via
#' [igraph::plot.igraph()], or dynamically via
#' [visNetwork::visIgraph()][visNetwork::visNetwork-igraph]
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
#' be displayed
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be displayed.
#' @param overlap_threshold Numeric value, between 0 and 1. Defines the threshold
#' to be used for removing edges in the enrichment map - edges below this value
#' will be excluded from the final graph. Defaults to 0.1.
#' @param scale_edges_width A numeric value, to define the scaling factor for the
#' edges between nodes. Defaults to 200 (works well chained to `visNetwork`
#' functions).
#' @param scale_nodes_size A numeric value, to define the scaling factor for the
#' node sizes. Defaults to 5 - works well chained to `visNetwork` functions.
#' @param color_by Character, specifying the column of `res_enrich` to be used
#' for coloring the plotted gene sets. Defaults to `gs_pvalue`.
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
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#' # annotation object
#' anno_df <- data.frame(
#'   gene_id = rownames(dds_macrophage),
#'   gene_name = mapIds(org.Hs.eg.db,
#'     keys = rownames(dds_macrophage),
#'     column = "SYMBOL",
#'     keytype = "ENSEMBL"
#'   ),
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
#'   res_de,
#'   anno_df,
#'   n_gs = 20
#' )
#'
#' em
#'
#' # could be viewed interactively with
#' # library("visNetwork")
#' # library("magrittr")
#' # em %>%
#' #   visIgraph() %>%
#' #   visOptions(highlightNearest = list(enabled = TRUE,
#' #                                      degree = 1,
#' #                                      hover = TRUE),
#' #             nodesIdSelection = TRUE)
enrichment_map_static <- function(res_enrich,
                           res_de,
                           annotation_obj,
                           gtl = NULL,
                           n_gs = 50,
                           gs_ids = NULL,
                           overlap_threshold = 0.1,
                           scale_edges_width = 200,
                           scale_nodes_size = 1,
                           color_by = "gs_pvalue",
                           cluster_fun = "cluster_markov") {

  cluster_fun <- match.arg(
    cluster_fun, c("cluster_markov", "cluster_louvain", "cluster_walktrap")
  )
  cluster_fun <- match.fun(cluster_fun)

  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }

  if (!color_by %in% colnames(res_enrich)) {
    stop(
      "Your res_enrich object does not contain the ",
      color_by,
      " column.\n",
      "Compute this first or select another column to use for the color."
    )
  }

  n_gs <- min(n_gs, nrow(res_enrich))

  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)], # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id] # the ones specified from the custom list
    )
  )
  
  overlap_matrix <- create_jaccard_matrix(res_enrich,
    n_gs = n_gs,
    gs_ids = gs_ids,
    return_sym = FALSE
  )

  rownames(overlap_matrix) <- colnames(overlap_matrix) <- res_enrich[rownames(overlap_matrix), "gs_description"]

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
  emg <- graph_from_data_frame(omm[, c(1, 2)], directed = FALSE)

  # add the edges, delete the ones that are under the threshold
  E(emg)$width <- sqrt(omm$value * scale_edges_width)
  emg <- delete_edges(emg, E(emg)[omm$value < overlap_threshold])

  # Find communities for layout
  gs_communities <- cluster_fun(emg)
  res_enrich$gs_membership <- factor(gs_communities$membership)
  V(emg)$membership <- gs_communities$membership
  # V(emg)$color <- gs_communities$membership

  idx <- match(V(emg)$name, res_enrich$gs_description)

  gs_size <- res_enrich$gs_de_count[idx]

  V(emg)$size <- scale_nodes_size * sqrt(gs_size)
  V(emg)$original_size <- gs_size

  # questo non mi serve perché li coloro in base al membership? oppure mi serve
  col_var <- res_enrich[idx, color_by]

  # the palette changes if it is z_score VS pvalue
  if (all(col_var <= 1) & all(col_var > 0)) { # likely p-values...
    col_var <- -log10(col_var)
    mypal <- colorRampPalette(RColorBrewer::brewer.pal(name = "YlOrRd", 9))(50)

    V(emg)$color <- mosdef::map_to_color(col_var, mypal, symmetric = FALSE, 
                                         limits = range(na.omit(col_var)))
   
    V(emg)$color[is.na(V(emg)$color)] <- "lightgrey"
  
  } else {
    # e.g. using z_score or aggregated value
    if (prod(range(na.omit(col_var))) >= 0) {
      # gradient palette
      mypal <- colorRampPalette(RColorBrewer::brewer.pal(name = "Oranges", 9))(50)

      V(emg)$color <- mosdef::map_to_color(col_var, mypal, symmetric = FALSE, 
                                           limits = range(na.omit(col_var)))
      
      V(emg)$color[is.na(V(emg)$color)] <- "lightgrey"
    } else {
      # divergent palette to be used
      mypal <- rev(colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50))
      V(emg)$color <- mosdef::map_to_color(col_var, mypal, symmetric = TRUE, 
                                           limits = range(na.omit(col_var)))
      V(emg)$color[is.na(V(emg)$color)] <- "lightgrey"
    }
  }

  V(emg)$color.border <- "black"

  # additional specification of edge colors
  E(emg)$color <- "lightgrey"

  # re-sorting the vertices alphabetically
  rank_gs <- rank(V(emg)$name)
  emg <- permute(emg, rank_gs)


  ground_truth = igraph::make_clusters(emg,V(emg)$membership)
  # plot(ground_truth,emg)
  # we need to change the layout
  emg_layout <- emg
  # emg_original <- emg

  browser()

  E(emg_layout)$weight <- apply(igraph::as_edgelist(emg_layout), 1, function(row) {
    weight.community(as.character(row), igraph::membership(gs_communities), 20, 1)
  })

  emg_layout$layout=igraph::layout_with_fr(emg_layout,weights=E(emg_layout)$weight)
  emg$layout=igraph::layout_with_fr(emg_layout,weights=E(emg_layout)$weight)
  
  cluster_centers <- sapply(igraph::groups(gs_communities), function(nodes) {
    node_indices <- match(nodes, V(emg_layout)$name) 
    centroid <- colMeans(igraph::layout_with_fr(emg_layout)[node_indices, , drop = FALSE])
    return(centroid)
  })

  cluster_labels <- add_cluster_names(gs_communities)
  V(emg)$label <- NA

  plot(emg, 
    mark.groups = igraph::communities(gs_communities)) 
  # Add word cloud labels at cluster centroids
  text(
    cluster_centers[, 1], 
    cluster_centers[, 2], 
    abels = cluster_labels, 
    col = "black", cex = 1.2, font = 2)

  # # add backbone links as edge attribute
  # plot(emg_layout)
  # plot(emg_layout, vertex.color = V(emg)$color, vertex.size = V(emg)$size, edge.width = E(emg)$width, edge.color = E(emg)$color)
  # plot(gs_communities, emg, 
  #   vertex.color = V(emg)$color, 
  #   vertex.size = V(emg)$size, 
  #   edge.width = E(emg)$width, 
  #   edge.color = E(emg)$color, 
  #   mark.groups = igraph::communities(gs_communities))
 



  # plot(emg_original)
  #  # E(g)$col <- FALSE
  # # E(g)$col[bb$backbone] <- TRUE
  # g_tbl <- tidygraph::as_tbl_graph(emg)
  # g_tbl <- g_tbl %>% mutate(cluster = igraph::membership(gs_communities))
  # cluster_centers$cluster <- seq_len(nrow(cluster_centers))  # Add cluster IDs

  # ggraph::ggraph(g_tbl, layout = "fr") +
  #   ggraph::geom_edge_link(aes(alpha = 0.5), color = "gray") +  # Edges
  #   ggraph::geom_node_point(aes(color = as.factor(cluster)), size = 6) +  # Nodes colored by cluster
  #   ggplot2::geom_text(aes(x = X1, y = X2, label = paste0("Cluster ", cluster)), 
  #             data = cluster_centers, size = 6, fontface = "bold", color = "black") +  # Cluster labels
  #   ggplot2::theme_void() + 
  #   ggplot2::theme(legend.position = "none")

  return(emg)
}

weight.community <- function(row,membership,weigth.within,weight.between){
  if(as.numeric(membership[which(names(membership)==row[1])]) == 
    as.numeric(membership[which(names(membership)==row[2])])){
    weight=weigth.within
  }else{
    weight=weight.between
  }
  return(weight)
}


add_cluster_names <- function(gs_communities, n_words = 3) {

    # Extract node names for each community
  community_texts <- lapply(groups(gs_communities), function(nodes) {
    sentences <- V(g)$name[match(nodes, V(g)$name)]  # Ensure correct matching
    paste(sentences, collapse = " ")  # Combine all node names into one string per cluster
  })

  # Tokenize words and count frequency
  word_counts <- lapply(community_texts, function(text) {
    data.frame(word = unlist(str_split(text, "\\s+"))) %>% 
      count(word, sort = TRUE) %>%
      filter(!word %in% stop_words$word)  # Remove common stop words
  })

  # Generate cluster labels (top 3 words)
  cluster_labels <- sapply(word_counts, function(df) {
    paste(head(df$word, n_words), collapse = " ")  # Take the top 3 words
  })

}