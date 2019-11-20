gs_dendro <- function(res_enrich,
                      n_gs = nrow(res_enrich),
                      gs_ids = NULL,
                      gs_dist_type = "kappa", # alternatives
                      clust_method = "ward.D2",
                      color_leaves_by = "z_score"
) {
  n_gs <- min(n_gs, nrow(res_enrich))

  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)],  # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id]  # the ones specified from the custom list
    )
  )

  if (gs_dist_type == "kappa") {
    dmat <- create_kappa_matrix(res_enrich, n_gs, gs_ids)
  } else if (gs_dist_type == "jaccard") {
    dmat <- create_jaccard_matrix(res_enrich, n_gs, gs_ids,return_sym = TRUE)
  } else {
    dmat <- GeneTonic:::create_semsim_matrix(res_enrich, semsim_data = semsim_data, n_gs, gs_ids)
  }

  rownames(dmat) <- colnames(dmat) <- res_enrich[gs_to_use, "gs_description"]

  my_dend <- as.dist(1 - dmat) %>%
    hclust(method = clust_method)

  dend_idx <- my_dend$order # keep the sorted index vector for the leaves

  # setup color
  mypal <- rev(scales::alpha(colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 1))
  col_var <- res_enrich[gs_to_use, "z_score"]
  leaves_col <- map2color(col_var, mypal, limits = range(col_var))[dend_idx]

  # setup size
  size_var <- -log10(as.numeric(res_enrich[gs_to_use, "gs_pvalue"]))[dend_idx]
  leaves_size <- 2 * (size_var - min(size_var)) / (max(size_var) - min(size_var) + 1e-10) + 0.2

  # leaves: colored by Z score
  # size: colored by squared size?
  # branches: indicating the treecut clusters

  my.clusters <- unname(dynamicTreeCut::cutreeDynamic(my_dend,
                                                      distM = as.matrix(dmat),
                                                      minClusterSize = 4,
                                                      verbose = 0))

  # or use
  clust_pal <- RColorBrewer::brewer.pal(max(my.clusters),"Set1")
  clust_cols <- (clust_pal[my.clusters])[dend_idx]


  # TODO: split up and add if clauses?
  my_dend %>%
    as.dendrogram(hang = -1) %>%
    branches_attr_by_clusters(my.clusters[dend_idx], values = clust_pal) %>%
    set("leaves_pch", 19) %>%             # TODO: can I set the border to black? (and leave the fill to black?)
    set("leaves_cex", leaves_size) %>%    # or to use gs_size?
    set("leaves_col", leaves_col) %>%
    # set("labels_colors", leaves_col) %>%
    plot(horiz = TRUE)
    # as.ggdend() %>% ggplot() %>% plotly::ggplotly()

  # TODO: or return the dend object and plot it out? this would enable ggdend to step in
}
