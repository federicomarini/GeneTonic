#' Creates the plot of the enrichment map graph
#'
#' Generates a graph plot for the enrichment map, as created within `GeneTonic`
#' via `enrichment_map()`.
#'
#' @param emg An `igraph` object, ideally generated via `enrichment_map()` 
#' @param cluster_fun Character string, containing the name of the function
#' to be used for clustering the graph of the geneset similarities. Defaults to
#' "cluster_markov".
#' @param scale_edges_width A numeric value, to define the scaling factor for the
#' edges between nodes. Defaults to 200 (works well chained to `visNetwork`
#' functions).
#' @param scale_nodes_size A numeric value, to define the scaling factor for the
#' within GeneTonic
#' @return A `ggraph` object with the static representation of the enrichment 
#' map graph.
#' 
#' @importFrom ggraph ggraph
#' @importFrom ggforce geom_mark_hull
#' @importFrom igraph edge_attr_names vertex_attr_names cluster_louvain
#' cluster_walktrap sizes groups delete_vertices as_edgelist membership
#' layout_with_kk
#' @importFrom scales alpha
#' 
#' @seealso [enrichment_map()] is used to generate the input igraph object.
#' Also, [GeneTonic()] embeds an interactive visualization for the enrichment 
#' map (based on the `VisNetwork` package)
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
#'   n_gs = 100
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
  # edge_attrs <- igraph::edge_attr_names(emg)
  E(emg)$width_not_scaled <- E(emg)$width * 0.5
  
  # Moved here to highlight the difference between the emapplot static and not emapplot
  cluster_fun <- match.arg(
    cluster_fun, c("cluster_markov", "cluster_louvain", "cluster_walktrap")
  )
  cluster_fun <- match.fun(cluster_fun)

  # Find the clusters to update the layout, add it to the graph
  gs_communities <- cluster_fun(emg)
  V(emg)$membership <- gs_communities$membership
  
  message("GeneTonicInfo: found ", length(table(gs_communities$membership)),
          " clusters of genesets")

  # There is a bug in hull function: if there is a cluster of 1 element it crashes, i delete the clusters of only one element
  # They can be added to the graph without hull. (todo)
  ## TODO: this part will be removed once v0.5.0 of ggforce hits CRAN
  small_clusters <- which(igraph::sizes(gs_communities) <= 2)
  nodes_to_remove <- unlist(igraph::groups(gs_communities)[small_clusters])
  
  if (length(nodes_to_remove) > 0)
    message("GeneTonicInfo: removing ", sum(small_clusters > 0),
            " clusters of genesets with one or two members only",
            " (total: ", length(nodes_to_remove), " nodes removed)")
  
  emg <- igraph::delete_vertices(emg, nodes_to_remove)

  # We create a layout dummy object to use artifcially created 
  # emg <- emg
  E(emg)$weight <- apply(igraph::as_edgelist(emg), 1, function(row) {
    weight_community(as.character(row), igraph::membership(gs_communities), 5, 2)
  })

  # A dataframe with the annotation of each cluster is created 
  # annotation -> merge all gesets name from cluster, count words, delete connection words, take first 4
  cluster_labels <- add_cluster_names(emg, gs_communities)

  # The membership is added to the graph as factor
  V(emg)$membership <- as.factor(as.character(V(emg)$membership))
  V(emg)$cluster_label <- 
    as.factor(cluster_labels[match(V(emg)$membership, names(cluster_labels))])

  # print(names(igraph::vertex_attr(emg)))
  # print(names(igraph::edge_attr(emg)))

  mem_df <- data.frame(names = V(emg)$name,
                       membership = as.numeric(V(emg)$membership))

  message("GeneTonicInfo: plotting ", length(table(mem_df$membership)),
          " clusters of genesets")


  #### Layout experiments
  # browser()
  E(emg)$weigth_to_use_for_layout <- log(E(emg)$width) * E(emg)$weight

  # hist((E(emg)$weigth_to_use))
  # lay <- layout_by_cluster(emg, mem_df, layout = igraph::layout_with_fr)
  lay <- layout_by_cluster(g = emg, 
                           mem_df = mem_df, 
                           layout_to_use = igraph::layout_with_kk)
  # lay <- layout_by_cluster(emg, mem_df, layout = igraph::layout_with_graphopt)
  # lay <- layout_by_cluster(emg, mem_df, layout = igraph::layout_with_graphopt)
  # lay <- igraph::layout_with_fr(emg) # ,weights=E(emg)$weight)
  # lay <- igraph::layout_with_graphopt(emg)
  # lay <- igraph::layout_with_kk(emg)
  # lay <- igraph::component_wise(lay)
  # lay <- igraph::layout_(emg, igraph::layout_with_kk(), igraph::component_wise()) # Source? https://igraph.org/r/html/1.3.0/layout_.html
  
  
  #### Artifical labels for the clusters
  # cluster_centers <- data.frame(t(sapply(igraph::groups(gs_communities), function(nodes) {
  #   node_indices <- match(nodes, V(emg)$name) 
  #   centroid <- colMeans(lay[node_indices, , drop = FALSE])
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

  gp <- 
    ggraph::ggraph(emg,
       layout = "manual",
       x = lay[, 1],
       y = lay[, 2]) +
    
    ## edges first, so they don't cover anything
    ggraph::geom_edge_link0(
      aes(edge_width = .data$width_not_scaled), 
      ## adding some transparency here
      edge_colour = scales::alpha("lightgrey", 0.7)) +
    
    ## hull on top, so that the nodes still are in the "native color"
    ggforce::geom_mark_hull(
      aes(.data$x, .data$y,
          fill= .data$cluster_label,
          ## why not having the border of the hull too "colored in sync"
          color = .data$cluster_label, 
          label = .data$cluster_label
          ),
      label.fill = scales::alpha("white", 0.1),
      concavity = 1000,
      expand = unit(3, "mm"),
      alpha = 0.15
    ) +
    
    ## handling the individual nodes at the end
    ggraph::geom_node_point(aes(size = .data$original_size, fill = I(V(emg)$color)), shape = 21, color = "black") + # I(V(emg_for_gggraph)$ # Adjust border thickness) +  # Use `I()` to prevent scaling
    # ggplot2::scale_fill_identity() +
    ggplot2::scale_size_continuous(range = c(7, 25)) +  # Adjust min and max sizes  
    
    ## background as clean as possible
    ggraph::theme_graph() + 
    ## possibly to be done optional?
    ggplot2::theme(legend.position = "none")
  
    # ggrepel::geom_label_repel(data = cluster_centers, 
    #                     aes(x = X1, y = X2, label = cluster_labels),
    #                     size = 2, 
    #                     color = "black", fill = alpha("white", .15), 
    #                     label.size = 0.5,  # Border thickness
    #                     label.padding = unit(0.2, "lines"))
      # ggplot2::geom_text(data = cluster_centers, 
      #       aes(x = X1, y = X2, label = cluster_labels),
      #       size = 3, 
      #       color = "black", 
      #       fontface = "bold")
  
  return(gp)
}

# To create the artificial weights
weight_community <- function(row, membership, weigth_within, weight_between) {
  if(as.numeric(membership[which(names(membership) == row[1])]) == 
     as.numeric(membership[which(names(membership) == row[2])])){
    weight <- weigth_within
  }else{
    weight <- weight_between
  }
  
  return(weight)
}

# To create the names of the clusters
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

# For this also component_wise from igraph seems to be useful? but i could not make it work.https://igraph.org/r/html/1.3.0/component_wise.html
#' Layout a graph, by cluster/community
#'
#' @param g An `igraph` graph object
#' @param mem_df A data frame containing the information on the communities
#' detected in the graph. Has to contain a column called "membership", and can 
#' be commonly derived by a `communities` object (e.g. output from 
#' `cluster_louvain()` or similar)
#' @param layout_to_use Function to use for laying out each subgraph. Defaults
#' to `igraph::layout_with_kk()`
#' 
#' @importFrom igraph induced_subgraph merge_coords layout_with_kk
#' disjoint_union
#'
#' @returns A matrix containing the coordinates, to be passed and set manually
#' for plotting
#' 
#' @details
#' This function is a re-implementation of the `layoutByCluster` function within
#' the BioNAR package. An additional parameter to easily customize the layout 
#' within each subgraph has been added, for higher flexibility.
#' 
#' 
#' @noRd
#'
#' @examples
#' # TODO - but it is an internal function per se...
layout_by_cluster <- function(g, 
                              mem_df, 
                              layout_to_use = igraph::layout_with_kk) {
    
    subgraph_list <- lapply(
      names(table(mem_df$membership)),
      function(cluster_id) {
        mem <- mem_df$membership
        idx <- which(mem == cluster_id)
        sg <- igraph::induced_subgraph(g, V(g)[idx])
        return(sg)
      }
    )
                      
    layout_list <- lapply(subgraph_list, function(sg) {
      if ((identical(layout_to_use, igraph::layout_with_fr) || identical(layout_to_use, igraph::layout_with_kk)) && "weight_to_use_for_layout" %in% E(sg)) {
      layout_to_use(sg, weights = E(sg)$weigth_to_use_for_layout)
      } else {
      layout_to_use(sg)
      }
    })

    layout_merged <- igraph::merge_coords(subgraph_list, layout_list)
    unified_graph <- igraph::disjoint_union(subgraph_list)
    idx <- match(V(g)$name, V(unified_graph)$name)
    lay <- layout_merged[idx, ]
    
    return(lay)
}

  # emg$layout <- igraph::layout_with_fr(emg,weights=E(emg)$weight)
  # emg$layout <- igraph::layout_with_fr(emg,weights=E(emg)$weight)

  # layout <- igraph::layout_with_fr(emg,weights=E(emg)$weight)


  # V(emg_for_gggraph)$membership <- as.factor(as.character(V(emg_for_gggraph)$membership))

  # cluster_labels <- add_cluster_names(emg ,gs_communities)
  # V(emg_for_gggraph)$cluster_label <- as.factor(cluster_labels[match(V(emg)$membership, names(cluster_labels))])
  # # cluster_labels <- names(igraph::groups(gs_communities))
  # V(emg_for_gggraph)$label <- NA
  # cluster_centers <- data.frame(t(sapply(igraph::groups(gs_communities), function(nodes) {
  #   node_indices <- match(nodes, V(emg)$name) 
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
  # lay <-BioNAR::layout_by_cluster(emg, mem_df, layout_to_use = igraph::layout_with_kk)

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
  # plot(emg)
  # plot(emg, vertex.color = V(emg)$color, vertex.size = V(emg)$size, edge.width = E(emg)$width, edge.color = E(emg)$color)
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
