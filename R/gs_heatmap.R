#' Plot a heatmap of the gene signature on the data
#'
#' Plot a heatmap for the selected gene signature on the provided data, with the possibility to compactly display also DE only genes
#'
#' @param se A `SummarizedExperiment` object, or an object derived from this class,
#' such as a `DESeqTransform` object (variance stabilized transformed data, or
#' regularized logarithm transformed), in where the transformation has been applied
#' to make the data more homoscedastic and thus a better fit for visualization.
#' @param res_de A `DESeqResults` object.
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#' @param geneset_id Character specifying the gene set identifier to be plotted
#' @param genelist A vector of character strings, specifying the identifiers
#' contained in the row names of the `se` input object.
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param FDR Numeric value, specifying the significance level for thresholding
#' adjusted p-values. Defaults to 0.05.
#' @param de_only Logical, whether to include only differentially expressed genes
#' in the plot
#' @param cluster_rows Logical, determining if rows should be clustered, as
#' specified by [ComplexHeatmap::Heatmap()]
#' @param cluster_columns Logical, determining if columns should be clustered, as
#' specified by [ComplexHeatmap::Heatmap()]
#' @param center_mean Logical, whether to perform mean centering on the row-wise
#' @param scale_row Logical, whether to standardize by row the expression values
#' @param winsorize_threshold Numeric value, to be applied as value to winsorize 
#' the extreme values of the heatmap. Should be a positive number. Defaults to 
#' NULL, which corresponds to not applying any winsorization. Suggested values: 
#' enter 2 or 3 if using row-standardized values (`scale_row` is TRUE), or visually
#' inspect the range of the values if using simply mean centered values.
#' @param anno_col_info A character vector of names in `colData(dds)` to use for
#' decorating the heatmap as annotation.
#' @param plot_title Character string, to specify the title of the plot,
#' displayed over the heatmap. If left to `NULL` as by default, it tries to use
#' the information on the geneset identifier provided
#' @param ...  Additional arguments passed to other methods, e.g. in the call to
#' [ComplexHeatmap::Heatmap()]
#'
#' @return A plot returned by the [ComplexHeatmap::Heatmap()] function
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
#' vst_macrophage <- vst(dds_macrophage)
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
#' gs_heatmap(vst_macrophage,
#'   res_de,
#'   res_enrich,
#'   anno_df,
#'   geneset_id = res_enrich$gs_id[1],
#'   cluster_columns = TRUE,
#'   anno_col_info = "condition"
#' )
gs_heatmap <- function(se,
                       res_de,
                       res_enrich,
                       annotation_obj = NULL,
                       gtl = NULL,
                       geneset_id = NULL,
                       genelist = NULL,
                       FDR = 0.05,
                       de_only = FALSE,
                       cluster_rows = TRUE,
                       cluster_columns = FALSE,
                       center_mean = TRUE,
                       scale_row = FALSE,
                       winsorize_threshold = NULL,
                       anno_col_info = NULL,
                       plot_title = NULL,
                       ...) {
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }
  
  if (!is.null(winsorize_threshold)) {
    stopifnot(is.numeric(winsorize_threshold))
    stopifnot(winsorize_threshold >= 0)
  }

  # check that the data would ideally be a DST, so that it is not the counts/normalized?
  mydata <- assay(se)

  # if(geneset in the results)
  #   pick the genes from there
  # else
  #   option to override the geneset by providing a list
  #
  # idea: multiselect with gene names - but in the UI
  # internal matching to the IDs (in this function we use the ids already)

  # rownames(res_enrich) <- res_enrich[["gs_id"]]
  if (!is.null(geneset_id)) {
    if (geneset_id %in% res_enrich[["gs_id"]]) {
      thisset_name <- res_enrich[geneset_id, "gs_description"]
      thisset_members <- unlist(strsplit(res_enrich[geneset_id, "gs_genes"], ","))
      thisset_members_ids <- annotation_obj$gene_id[match(thisset_members, annotation_obj$gene_name)]
    }
  } else {
    # overridable via a list
    if (!all(genelist %in% rownames(se))) {
      not_there <- genelist[!(genelist %in% rownames(se))]
      warning(
        "Some of the provided gene ids were not found in the SummarizedExperiment",
        "\nNot found: ", not_there
      )
    }
    thisset_members_ids <- intersect(rownames(se), genelist)
    thisset_name <- "Custom list"
  }

  sig_to_keep <- (thisset_members_ids %in% rownames(se)) #
  thisset_members_ids_available <- thisset_members_ids[sig_to_keep]

  mydata_sig <- mydata[thisset_members_ids_available, ]

  # to avoid problems later, remove the ones non-expressed and with variance = 0
  to_remove <- apply(mydata_sig, 1, var) == 0
  mydata_sig <- mydata_sig[!to_remove, ]

  hm_name <- "Expression \nvalues"

  if (center_mean) {
    mydata_sig <- mydata_sig - rowMeans(mydata_sig)
    hm_name <- "Expression \nvalues"
  }

  if (scale_row) {
    mydata_sig <- t(scale(t(mydata_sig)))
    hm_name <- "Z-scores \nExpression \nvalues"
  }

  if (de_only) {
    de_res <- deseqresult2df(res_de, FDR)
    de_genes <- de_res$id
    de_to_keep <- rownames(mydata_sig) %in% de_genes
    mydata_sig <- mydata_sig[de_to_keep, , drop = FALSE]
  }

  extreme_value <- max(abs(range(mydata_sig)))
  if (!is.null(winsorize_threshold)) {
    # do the winsoring
    mydata_sig[mydata_sig < -winsorize_threshold] <- -winsorize_threshold
    mydata_sig[mydata_sig > winsorize_threshold] <- winsorize_threshold
  } 
  # dim(mydata_sig)

  if (is.null(plot_title)) {
    title <- paste0("Signature heatmap - ", thisset_name, " - ", geneset_id)
  } else {
    title <- plot_title
  }

  ### anno_col_info <- anno_col_info[anno_col_info %in% colnames(colData(se))]
  ### sample_decoration <- as.data.frame(colData(se))[, anno_col_info, drop = FALSE]

  # could there be a way to make this programmatically & clever?

  ## if only one column: SO
  # anno_col_vals <- colData(se)[,anno_col_info,drop = TRUE]
  #
  # ha_cols <- list(
  #   Annotation = structure(
  #     brewer.pal(length(unique(anno_col_vals)), "Set1"),
  #     names = unique(as.character(anno_col_vals))
  #   )
  # )
  # deco_ha <- HeatmapAnnotation(
  #   name = "eheh",
  #   Annotation = anno_col_vals,
  #   col = ha_cols
  # )
  ### deco_ha <- HeatmapAnnotation(df = sample_decoration)
  # if(returnData) {
  #
  # }

  # pheatmap(mydata_sig,
  #          # annotation_col = anno_colData,
  #          cluster_rows = cluster_rows, cluster_cols = cluster_cols,
  #          scale = ifelse(scale_row, "row", "none"), main = title,
  #          labels_row = annotation_obj[rownames(mydata_sig), ]$gene_name,
  #          annotation_col = sample_decoration)

  if (is.null(anno_col_info)) {
    ch <- ComplexHeatmap::Heatmap(
      matrix = mydata_sig,
      column_title = title,
      name = hm_name,
      col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
      rect_gp = gpar(col = "white", lwd = 0.5),
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      row_labels = annotation_obj$gene_name[match(rownames(mydata_sig), annotation_obj$gene_id)],
      ...
    )
  } else {
    anno_col_info <- anno_col_info[anno_col_info %in% colnames(colData(se))]
    sample_decoration <- as.data.frame(colData(se))[, anno_col_info, drop = FALSE]

    # handling the color with user defined values
    col_list <- vector(mode = "list", length = ncol(sample_decoration))

    for (i in seq_len(ncol(sample_decoration))) {
      these_decos <- sample_decoration[, i, drop = TRUE]

      if (is.numeric(these_decos)) {
        max_val <- max(these_decos, na.rm = TRUE)
        min_val <- min(these_decos, na.rm = TRUE)
        mid_val <- median(these_decos, na.rm = TRUE)

        if (max_val * min_val >= 0) {
          # sequential palette, all the same sign
          these_cols <- circlize::colorRamp2(
            c(min_val, max_val),
            c("grey95", "darkred")
          )
        } else {
          # need for a diverging palette
          these_cols <- circlize::colorRamp2(
            c(min_val, mid_val, max_val),
            c("blue", "grey95", "red")
          )
        }
      }

      if (is.character(these_decos) | is.factor(these_decos)) {
        if (is.character(these_decos)) {
          these_decos <- factor(these_decos)
        }
        these_cols <- colorspace::rainbow_hcl(n = length(levels(these_decos)))
        names(these_cols) <- levels(these_decos)
      }

      col_list[[i]] <- these_cols
    }

    names(col_list) <- colnames(sample_decoration)

    deco_ha <- HeatmapAnnotation(df = sample_decoration, col = col_list)

    ch <- ComplexHeatmap::Heatmap(
      matrix = mydata_sig,
      column_title = title,
      name = hm_name,
      col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
      rect_gp = gpar(col = "white", lwd = 0.5),
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      row_labels = annotation_obj$gene_name[match(rownames(mydata_sig), annotation_obj$gene_id)],
      top_annotation = deco_ha,
      ...
    )
  }
  draw(ch, merge_legend = TRUE)
}



#' Compute gene set scores
#'
#' Compute gene set scores for each sample, by transforming the gene-wise change
#' to a geneset-wise change
#'
#' @param se A `SummarizedExperiment` object, or an object derived from this class,
#' such as a `DESeqTransform` object (variance stabilized transformed data, or
#' regularized logarithm transformed), in where the transformation has been applied
#' to make the data more homoscedastic and thus a better fit for visualization.
#' @param res_de A `DESeqResults` object.
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param annotation_obj A `data.frame` object with the feature annotation
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#'
#' @return A matrix with the geneset Z scores, e.g. to be plotted with [gs_scoresheat()]
#'
#' @seealso [gs_scoresheat()] plots these scores
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
#' vst_macrophage <- vst(dds_macrophage)
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
#' scores_mat <- gs_scores(
#'   vst_macrophage,
#'   res_de,
#'   res_enrich[1:50, ],
#'   anno_df
#' )
gs_scores <- function(se,
                      res_de, # maybe won't be needed?
                      res_enrich,
                      annotation_obj = NULL,
                      gtl = NULL) {
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }

  mydata <- assay(se)
  # returns a matrix, rows = genesets, cols = samples
  # rownames(res_enrich) <- res_enrich[["gs_id"]]

  rowsd_se <- matrixStats::rowSds(mydata)
  rowavg_se <- rowMeans(mydata)
  mydata_z <- (mydata - rowavg_se) / rowsd_se

  gss_mat <- matrix(NA, nrow = nrow(res_enrich), ncol = ncol(se))
  rownames(gss_mat) <- paste0(
    res_enrich[["gs_description"]], "|",
    res_enrich[["gs_id"]]
  )
  # rownames(gss_mat) <- res_enrich[["gs_id"]]
  colnames(gss_mat) <- colnames(se)

  for (i in seq_len(nrow(res_enrich))) {
    thisset_members <- unlist(strsplit(res_enrich[i, "gs_genes"], ","))
    thisset_members_ids <- annotation_obj$gene_id[match(thisset_members, annotation_obj$gene_name)]

    thisset_members_ids <- thisset_members_ids[thisset_members_ids %in% rownames(se)]
    thisset_zs <- mydata_z[thisset_members_ids, , drop = FALSE]
    thisset_mean_zs <- colMeans(thisset_zs)

    gss_mat[i, ] <- thisset_mean_zs
  }

  return(gss_mat)
}

#' Plots a matrix of geneset scores
#'
#' Plots a matrix of geneset Z scores, across all samples
#'
#' @param mat A matrix, e.g. returned by the [gs_scores()] function
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be displayed.
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be displayed.
#' @param clustering_distance_rows Character, a distance measure used in
#' clustering rows
#' @param clustering_distance_cols Character, a distance measure used in
#' clustering columns
#' @param cluster_rows Logical, determining if rows should be clustered
#' @param cluster_cols Logical, determining if columns should be clustered
#'
#' @return A `ggplot` object
#'
#' @seealso [gs_scores()] computes the scores plotted by this function
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
#' vst_macrophage <- vst(dds_macrophage)
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
#' scores_mat <- gs_scores(
#'   vst_macrophage,
#'   res_de,
#'   res_enrich[1:30, ],
#'   anno_df
#' )
#' gs_scoresheat(scores_mat,
#'   n_gs = 30
#' )
gs_scoresheat <- function(mat,
                          n_gs = nrow(mat),
                          gs_ids = NULL,
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean",
                          cluster_rows = TRUE,
                          cluster_cols = TRUE) {
  n_gs <- min(n_gs, nrow(mat))
  current_ids <- unlist(lapply(
    strsplit(rownames(mat), split = "|", fixed = TRUE), function(arg) arg[[2]]
  ))

  gs_to_use <- unique(
    c(
      current_ids[seq_len(n_gs)], # the ones from the top
      gs_ids[gs_ids %in% current_ids] # the ones specified from the custom list
    )
  )

  gs_to_use <- match(gs_to_use, current_ids)
  mat <- mat[gs_to_use, ]

  d_rows <- dist(mat, method = clustering_distance_rows)
  d_cols <- dist(t(mat), method = clustering_distance_cols)

  row_tree <- hclust(d_rows)
  col_tree <- hclust(d_cols)
  score_matrix <- mat

  if (cluster_rows) {
    score_matrix <- score_matrix[row_tree$order, ]
  }

  if (cluster_cols) {
    score_matrix <- score_matrix[, col_tree$order]
  }

  labels_rows <- factor(rownames(score_matrix),
    levels = rev(rownames(score_matrix))
  )
  # to have the top ones on top, if not clustered
  labels_cols <- factor(colnames(score_matrix),
    levels = colnames(score_matrix)
  )

  score_df <- expand.grid(list(labels_rows, labels_cols),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  scores <- data.frame(as.vector(score_matrix))
  score_df <- cbind(score_df, scores)
  colnames(score_df) <- c("GeneSet", "Sample", "Value")

  p <- ggplot(score_df, aes_string(x = "Sample", y = "GeneSet")) +
    geom_tile(aes_string(fill = "Value"), color = "white") +
    # from brewer.pal(n = 11, name = 'RdYlBu') and brewer.pal(n = 11, name = 'YlGn') - or Spectal, 11?
    scale_fill_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(
        angle = 60,
        hjust = 1
      ),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 11)
    ) +
    labs(
      fill = "Geneset \nZ value"
    )
  return(p)
}
