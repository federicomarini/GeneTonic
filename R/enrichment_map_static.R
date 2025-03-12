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
#' @importFrom ggraph ggraph
#' @importFrom ggforce geom_mark_hull
#' @importFrom igraph edge_attr_names vertex_attr_names cluster_louvain
#' cluster_walktrap sizes groups delete_vertices as_edgelist membership
#' layout_with_kk
#' @importFrom scales alpha
#' @importFrom BioNAR layoutByCluster
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
#'   n_gs = 50
#' )
#' 
#' em2 <- enrichment_map(res_enrich,
#'   res_de,
#'   anno_df, color_by = "z_score",
#'   n_gs = 50
#' )
#'
#' em
#' 
#' ## TODO now works with...
#' gtl <- GeneTonicList(dds_macrophage, 
#'                      res_de = res_macrophage_IFNg_vs_naive, 
#'                      res_enrich = res_enrich, 
#'                      annotation_obj = anno_df)
#' ### enrichment_map_static(gtl = gtl, cluster_fun = "cluster_louvain")
#'                      
#' ## TODO: this should work with (for example)
#' plot_emap_static(em, cluster_fun = "cluster_markov")
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
plot_emap_static <- function(emg,
                             scale_edges_width = 5,
                             scale_nodes_size = 10,
                             cluster_fun = "cluster_markov") {

  stopifnot(is(emg, "igraph"))
  
  what_attrs_are_in <- vertex_attr_names(emg)
  
  
  ## TODO: We would ideally need to start from here, possibly "just transferring the"
  ## colors from the "interactive schemes"

  # The weights are not scaled in the same way (this can also be improved)
  V(emg)$weight <- V(emg)$original_size
  
  V(emg)$color <- V(emg)$color.highlight
  
  
  ## TODO: this would ideally require a different number of values?
  ## TODO: in the example, E(emg) is 199 edges, but omm$value is 
  #### E(emg)$width_not_scaled <- omm$value * scale_edges_width
  edge_attrs <- igraph::edge_attr_names(emg)
  
  E(emg)$width_not_scaled <- E(emg)$width * 0.5
  
  # Moved here to highlight the difference between the emapplot static and not emapplot
  cluster_fun <- match.arg(
    cluster_fun, c("cluster_markov", "cluster_louvain", "cluster_walktrap")
  )
  cluster_fun <- match.fun(cluster_fun)

  # Find the clusters to update the layout, add it to the graph
  gs_communities <- cluster_fun(emg)
  V(emg)$membership <- gs_communities$membership

  # There is a bug in hull function: if there is a cluster of 1 element it crashes, i delete the clusters of only one element
  # They can be added to the graph without hull. (todo)
  ## TODO: this part will be removed once v0.5.0 of ggforce hits CRAN
  small_clusters <- which(igraph::sizes(gs_communities) <= 2)
  nodes_to_remove <- unlist(igraph::groups(gs_communities)[small_clusters])
  emg <- igraph::delete_vertices(emg, nodes_to_remove)

  # We create a layout dummy object to use artifcially created 
  emg_layout <- emg
  E(emg_layout)$weight <- apply(igraph::as_edgelist(emg_layout), 1, function(row) {
    weight.community(as.character(row), igraph::membership(gs_communities), 3, 1)
  })

  # A dataframe with the annotation of each cluster is created 
  # annotation -> merge all gesets name from cluster, count words, delete connection words, take first 4
  cluster_labels <- add_cluster_names(emg, gs_communities)

  # The membership is added to the graph as factor
  V(emg)$membership <- as.factor(as.character(V(emg)$membership))
  V(emg)$cluster_label <- 
    as.factor(cluster_labels[match(V(emg)$membership, names(cluster_labels))])

  mem.df <- data.frame(names = V(emg)$name,membership = as.numeric(V(emg)$membership))
  lay <- BioNAR::layoutByCluster(emg_layout, mem.df, layout = igraph::layout_with_kk)

  message("GeneTonicInfo: found ", length(table(mem.df$membership)),
          " clusters of genesets")
  
  # lay <- igraph::layout_with_fr(emg_layout,weights=E(emg_layout)$weight)
  
  gp <- ggraph::ggraph(emg,
       layout = "manual",
       x = lay[, 1],
       y = lay[, 2]) +
    
    ## edges first, so they don't cover anything
    ggraph::geom_edge_link0(
      aes(edge_width = width_not_scaled), 
      ## adding some transparency here
      edge_colour = scales::alpha("lightgrey", 0.7)) +
    
    ## hull on top, so that the nodes still are in the "native color"
    ggforce::geom_mark_hull(
      aes(x, y,
          fill= cluster_label,
          ## why not having the border of the hull too "colored in sync"
          color = cluster_label,
          label = cluster_label),
      label.fill = scales::alpha("white", 0.1),
      concavity = 10,
      expand = unit(7, "mm"),
      alpha = 0.15
    ) +
    
    ## handling the individual nodes at the end
    ggraph::geom_node_point(aes(size = size, fill = I(V(emg)$color)), shape = 21, color = "black") + # I(V(emg_for_gggraph)$ # Adjust border thickness) +  # Use `I()` to prevent scaling
    # ggplot2::scale_fill_identity() +
    ggplot2::scale_size_continuous(range = c(7, 25)) +  # Adjust min and max sizes  
    
    ## background as clean as possible
    ggraph::theme_graph() + 
    ## possibly to be done optional?
    ggplot2::theme(legend.position = "none")
  
    # ggrepel::geom_label_repel(data = cluster_centers, 
    #                     aes(x = X1, y = X2, label = cluster_labels),
    #                     size = 5, fontface = "bold", 
    #                     color = "black", fill = alpha("white", .15), 
    #                     label.size = 0.5,  # Border thickness
    #                     label.padding = unit(0.2, "lines"))+
  
  
  return(gp)
}

weight.community <- function(row, membership, weigth.within, weight.between) {
  if(as.numeric(membership[which(names(membership) == row[1])]) == 
     as.numeric(membership[which(names(membership) == row[2])])){
    weight <- weigth.within
  }else{
    weight <- weight.between
  }
  
  return(weight)
}


add_cluster_names <- function(emg, gs_communities, n_words = 4) {

  stop_words <- c("the", "and", "is", "in", "to", "of", 
                  "it", "that", "on", "for", "with", 
                  "as", "this", "was", "but", "be", 
                  "by", "or", "not", "are", "at", "an")
  
  # Extract node names for each community
  community_texts <- lapply(igraph::groups(gs_communities), function(nodes) {
    sentences <- igraph::V(emg)$name[match(nodes, igraph::V(emg)$name)]  # Ensure correct matching
    paste(sentences, collapse = " ")  # Combine all node names into one string per cluster
  })
  
  # Tokenize words and count frequency
  word_counts <- lapply(community_texts, function(text) {
    # dplyr::tibble(word = as.character(unlist(strsplit(text, "\\s+")))) %>% 
    #       dplyr::count(word, sort = TRUE) %>%
    #       dplyr::filter(!word %in% tidytext::stop_words$word)  # Remove common stop words
    
    words <- unlist(strsplit(text, "\\s+"))  # Split text into words
    # words <- words[!(words %in% tidytext::stop_words$word)]  # Remove stop words
    words <- words[!(words %in% stop_words)]  # Remove stop words
    
    # Create frequency table
    word_freq <- table(words)
    word_freq <- sort(word_freq, decreasing = TRUE)  # Sort by frequency
    
    # Convert to data frame
    df <- data.frame(word = names(word_freq), freq = as.numeric(word_freq), stringsAsFactors = FALSE)
    
    return(df)
  })
  
  # Generate cluster labels (top 3 words)
  cluster_labels <- sapply(word_counts, function(df) {
    paste(head(df$word, n_words), collapse = "\n")  # Take the top 3 words
  })

  return(cluster_labels)
}

  # emg_layout$layout <- igraph::layout_with_fr(emg_layout,weights=E(emg_layout)$weight)
  # emg$layout <- igraph::layout_with_fr(emg_layout,weights=E(emg_layout)$weight)

  # layout <- igraph::layout_with_fr(emg_layout,weights=E(emg_layout)$weight)


  # V(emg_for_gggraph)$membership <- as.factor(as.character(V(emg_for_gggraph)$membership))

  # cluster_labels <- add_cluster_names(emg ,gs_communities)
  # V(emg_for_gggraph)$cluster_label <- as.factor(cluster_labels[match(V(emg)$membership, names(cluster_labels))])
  # # cluster_labels <- names(igraph::groups(gs_communities))
  # V(emg_for_gggraph)$label <- NA
  # cluster_centers <- data.frame(t(sapply(igraph::groups(gs_communities), function(nodes) {
  #   node_indices <- match(nodes, V(emg_layout)$name) 
  #   centroid <- colMeans(layout[node_indices, , drop = FALSE])
  #   return(centroid)
  # })))

  # cluster_annotation <- merge(cluster_centers, as.data.frame(cluster_labels), by.x = "row.names", by.y = "row.names")
  # igraph::vertex_attr_names(emg)
  # community_colors <- rainbow(length(igraph::groups(gs_communities)))
  # centroid_df <- data.frame(
  #   x = cluster_centers[1, ],
  #   y = cluster_centers[2, ],
  #   cluster_label = cluster_labels[match(colnames(cluster_centers), names(cluster_labels))]
  # )
  # Remove clusters with 2 or fewer nodes
  # lay <-BioNAR::layoutByCluster(emg, mem.df, layout = igraph::layout_with_kk)

  # ggplot2::geom_text(
  #   data = centroid_df,  # Add the centroids as a data source
  #   aes(label = cluster_label, x = x, y = y),
  #   fontface = "bold", size = 3, vjust = 1.5, hjust = 0.5  # Adjust label positioning
  # )

  # ggplot2::scale_color_brewer(palette = "Set1") +
  # ggplot2::scale_fill_brewer(palette = "Set1") +
  # ggraph::scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
  # ggraph::theme_graph() +
  # ggplot2::theme(legend.position = "none")

  # Here and under is the part that needs to work, convverting the object to gggraph and then plotting it!

  # Convert igraph  object to tidygraph object
  # g_tbl <- tidygraph::as_tbl_graph(emg)

  # # Add cluster membership as a node attribute
  # g_tbl <- g_tbl %>% mutate(cluster = igraph::membership(gs_communities))

  # # Create a gggraph plot
  # ggraph::ggraph(g_tbl, layout = "fr") +
  #   ggraph::geom_edge_link(aes(alpha = 0.5), color = "gray") +  # Edges
  #   ggraph::geom_node_point(aes(color = as.factor(cluster)), size = 6) +  # Nodes colored by cluster
  #   ggplot2::geom_text(aes(x = cluster_centers[1, ], y = cluster_centers[2, ], label = cluster_labels), 
  #   #           size = 6, fontface = "bold", color = "black") +  # Cluster labels
  #   ggplot2::theme_void() + 
  #   ggplot2::theme(legend.position = "none")

  # plot(emg, 
  #   mark.groups = igraph::communities(gs_communities)) 
  # # Add word cloud labels at cluster centroids
  # text(
  #   cluster_centers[, 1], 
  #   cluster_centers[, 2], 
  #   labels = cluster_labels, 
  #   col = "black", cex = 1.2, font = 2)
  # Add legend for community colors

  # legend("topright",  # Position of the legend
  #       legend = cluster_labels,  # Community labels
  #       fill = community_colors,  # Colors corresponding to the communities
  #       border = "black",  # Border color for the legend boxes
  #       bty = "n",  # No box around the legend
  #       title = "Communities",  # Legend title
  #       cex = 0.8)  # Adjust legend text size
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

  # community_graph <- igraph::groups(gs_communities)
