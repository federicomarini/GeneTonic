#' Plots a summary of enrichment results
#'
#' Plots a summary of enrichment results for one set
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be displayed
#' @param p_value_column Character string, specifying the column of `res_enrich`
#' where the p-value to be represented is specified. Defaults to `gs_pvalue`
#' (it could have other values, in case more than one p-value - or an adjusted
#' p-value - have been specified).
#' @param color_by Character, specifying the column of `res_enrich` to be used
#' for coloring the plotted gene sets. Defaults sensibly to `z_score`.
#' @param return_barchart Logical, whether to return a barchart (instead of the 
#' default dot-segment plot); defaults to FALSE.
#'
#' @return A `ggplot` object
#'
#' @seealso [gs_summary_overview_pair()], [gs_horizon()]
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
#'
#' # res_enrich object
#' data(res_enrich_macrophage, package = "GeneTonic")
#' res_enrich <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
#' res_enrich <- get_aggrscores(res_enrich, res_de, anno_df)
#'
#' gs_summary_overview(res_enrich)
#' 
#' # if desired, it can also be shown as a barplot
#' gs_summary_overview(res_enrich, 30, return_barchart = TRUE)
gs_summary_overview <- function(res_enrich,
                                n_gs = 20,
                                p_value_column = "gs_pvalue",
                                color_by = "z_score",
                                return_barchart = FALSE
                                # , size_by = "gs_de_count"
                                ) {
  if (!is.null(color_by)) {
    if (!(color_by %in% colnames(res_enrich))) {
      stop("Your res_enrich object does not contain the ",
           color_by,
           " column.\n",
           "Compute this first or select another column to use for the color.")
    }
  }

  re <- res_enrich
  re$logp10 <- -log10(res_enrich[[p_value_column]])
  re <- re[seq_len(n_gs), ]

  re_sorted <- re %>%
    arrange(.data$logp10) %>%
    mutate(gs_description = factor(.data$gs_description, .data$gs_description))
  
  if (return_barchart) {
    p <- ggplot(re_sorted, (aes_string(x = "gs_description", y = "logp10"))) 
    
    if (is.null(color_by)) {
      p <- p + geom_col()
    } else {
      p <- p + geom_col(aes(fill = .data[[color_by]])) + 
        scale_fill_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026")
    }
    
    p <- p +
      coord_flip() + 
      labs(x = "Gene set description",
           y = "-log10 p-value",
           col = color_by) +
      theme_minimal()
  } else {
    p <- ggplot(re_sorted, (aes_string(x = "gs_description", y = "logp10"))) +
      geom_segment(aes_string(x = "gs_description", 
                              xend = "gs_description", 
                              y = 0, 
                              yend = "logp10"), color = "grey")
    
    if (is.null(color_by)) {
      p <- p + geom_point(size = 4)
    } else {
      p <- p + geom_point(aes(col = .data[[color_by]]), size = 4) +
        scale_color_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026")
    }
    
    p <- p +
      coord_flip() +
      labs(x = "Gene set description",
           y = "-log10 p-value",
           col = color_by) +
      theme_minimal()
  }
  return(p)
}

#' Plots a summary of enrichment results
#'
#' Plots a summary of enrichment results - for two sets of results
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param res_enrich2 As `res_enrich`, the result of functional enrichment analysis,
#' in a scenario/contrast different than the first set.
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be displayed
#' @param p_value_column Character string, specifying the column of `res_enrich`
#' where the p-value to be represented is specified. Defaults to `gs_pvalue`
#' (it could have other values, in case more than one p-value - or an adjusted
#' p-value - have been specified).
#' @param color_by Character, specifying the column of `res_enrich` to be used
#' for coloring the plotted gene sets. Defaults sensibly to `z_score`.
#' @param alpha_set2 Numeric value, between 0 and 1, which specified the alpha
#' transparency used for plotting the points for gene set 2.
#'
#' @return A `ggplot` object
#'
#' @seealso [gs_summary_overview()], [gs_horizon()]
#'
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
#' res_enrich2 <- res_enrich[1:42, ]
#' set.seed(42)
#' shuffled_ones <- sample(seq_len(42)) # to generate permuted p-values
#' res_enrich2$gs_pvalue <- res_enrich2$gs_pvalue[shuffled_ones]
#' res_enrich2$z_score <- res_enrich2$z_score[shuffled_ones]
#' res_enrich2$aggr_score <- res_enrich2$aggr_score[shuffled_ones]
#' # ideally, I would also permute the z scores and aggregated scores
#' gs_summary_overview_pair(res_enrich = res_enrich,
#'                          res_enrich2 = res_enrich2)
gs_summary_overview_pair <- function(res_enrich,
                                     res_enrich2,
                                     n_gs = 20,
                                     p_value_column = "gs_pvalue",
                                     color_by = "z_score",
                                     alpha_set2 = 1) {
  if (!(color_by %in% colnames(res_enrich))) {
    stop("Your res_enrich object does not contain the ",
         color_by,
         " column.\n",
         "Compute this first or select another column to use for the color.")
  }
  # same for set2
  if (!(color_by %in% colnames(res_enrich2))) {
    stop("Your res_enrich object does not contain the ",
         color_by,
         " column.\n",
         "Compute this first or select another column to use for the color.")
  }

  gs_set1 <- res_enrich$gs_id
  gs_set2 <- res_enrich2$gs_id
  gs_common <- intersect(gs_set1, gs_set2)

  if (length(gs_common) == 0) {
    stop("No gene sets have been found in common to the two enrichment results")
  }

  # restrict to the top common n_gs
  gs_common <- gs_common[seq_len(min(n_gs, length(gs_common)))]

  common_re1 <- res_enrich[gs_common, ]
  common_re2 <- res_enrich2[gs_common, ]

  common_re1$logp10 <- -log10(common_re1[[p_value_column]])
  common_re2$logp10 <- -log10(common_re2[[p_value_column]])

  re_both <- common_re1
  re_both[["logp10_2"]] <- common_re2$logp10
  re_both[[color_by]] <- common_re1[[color_by]]
  re_both[[paste0(color_by, "_2")]] <- common_re2[[color_by]]

  re_both_sorted <- re_both %>%
    arrange(.data$logp10) %>%
    mutate(gs_description = factor(.data$gs_description, .data$gs_description))

  p <- ggplot(re_both_sorted, aes_string(x = "gs_description", y = "logp10")) +
    geom_segment(aes_string(x = "gs_description", xend = "gs_description", y = "logp10_2", yend = "logp10"), color = "grey") +
    geom_point(aes(fill = .data[[color_by]]), size = 4, pch = 21) +
    geom_point(aes_string(y = "logp10_2", col = paste0(color_by, "_2")),
               size = 4, alpha = alpha_set2) +
    scale_color_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026", name = paste0(color_by, " set 2")) +
    scale_fill_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026", name = paste0(color_by, " set 1")) +
    coord_flip() +
    labs(x = "Gene set description",
         y = "-log10 p-value",
         col = color_by) +
    ylim(0, NA) +
    theme_minimal()

  return(p)
}


#' Plots a summary of enrichment results
#'
#' Plots a summary of enrichment results - horizon plot to compare one or more
#' sets of results
#'
#' @details It makes sense to have the results in `res_enrich` sorted by
#' increasing `gs_pvalue`, to make sure the top results are first sorted by the
#' significance (when selecting the common gene sets across the `res_enrich`
#' elements provided in `compared_res_enrich_list`)
#'
#' The gene sets included are a subset of the ones in common to all different
#' scenarios included in `res_enrich` and the elements of `compared_res_enrich_list`.
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param compared_res_enrich_list A named list, where each element is a `data.frame`
#' formatted like the standard `res_enrich` objects used by `GeneTonic`. The
#' names of the list are the names of the scenarios.
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be displayed
#' @param p_value_column Character string, specifying the column of `res_enrich`
#' where the p-value to be represented is specified. Defaults to `gs_pvalue`
#' (it could have other values, in case more than one p-value - or an adjusted
#' p-value - have been specified).
#' @param color_by Character, specifying the column of `res_enrich` to be used
#' for coloring the plotted gene sets. Defaults sensibly to `z_score`.
#' @param ref_name Character, defining the name of the scenario to compare
#' against (the one in `res_enrich`) - defaults to "ref_scenario".
#' @param sort_by Character string, either "clustered", or "first_set". This
#' controls the sorting order of the included terms in the final plot.
#' "clustered" presents the terms grouped by the scenario where they assume the
#' highest values. "first_set" sorts the terms by the significance value in the
#' reference scenario.
#'
#' @return A `ggplot` object
#'
#' @seealso [gs_summary_overview()], [gs_summary_overview_pair()]
#'
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
#' res_enrich2 <- res_enrich[1:42, ]
#' res_enrich3 <- res_enrich[1:42, ]
#' res_enrich4 <- res_enrich[1:42, ]
#'
#' set.seed(2*42)
#' shuffled_ones_2 <- sample(seq_len(42)) # to generate permuted p-values
#' res_enrich2$gs_pvalue <- res_enrich2$gs_pvalue[shuffled_ones_2]
#' res_enrich2$z_score <- res_enrich2$z_score[shuffled_ones_2]
#' res_enrich2$aggr_score <- res_enrich2$aggr_score[shuffled_ones_2]
#'
#' set.seed(3*42)
#' shuffled_ones_3 <- sample(seq_len(42)) # to generate permuted p-values
#' res_enrich3$gs_pvalue <- res_enrich3$gs_pvalue[shuffled_ones_3]
#' res_enrich3$z_score <- res_enrich3$z_score[shuffled_ones_3]
#' res_enrich3$aggr_score <- res_enrich3$aggr_score[shuffled_ones_3]
#'
#' set.seed(4*42)
#' shuffled_ones_4 <- sample(seq_len(42)) # to generate permuted p-values
#' res_enrich4$gs_pvalue <- res_enrich4$gs_pvalue[shuffled_ones_4]
#' res_enrich4$z_score <- res_enrich4$z_score[shuffled_ones_4]
#' res_enrich4$aggr_score <- res_enrich4$aggr_score[shuffled_ones_4]
#'
#' compa_list <- list(
#'   scenario2 = res_enrich2,
#'   scenario3 = res_enrich3,
#'   scenario4 = res_enrich4
#' )
#'
#' gs_horizon(res_enrich,
#'            compared_res_enrich_list = compa_list,
#'            n_gs = 50,
#'            sort_by = "clustered")
#' gs_horizon(res_enrich,
#'            compared_res_enrich_list = compa_list,
#'            n_gs = 20,
#'            sort_by = "first_set")
gs_horizon <- function(res_enrich,
                       compared_res_enrich_list,
                       n_gs = 20,
                       p_value_column = "gs_pvalue",
                       color_by = "z_score",
                       ref_name = "ref_scenario",
                       sort_by = c("clustered", "first_set")) {
  if (!(color_by %in% colnames(res_enrich))) {
    stop("Your res_enrich object does not contain the ",
         color_by,
         " column.\n",
         "Compute this first or select another column to use for the color.")
  }

  if (!n_gs > 0) {
    stop("Please select a value for `n_gs` greater than 0")
  }

  if (is.null(names(compared_res_enrich_list))) {
    message("You provided a list for comparison without specifying names, adding some defaults")
    names(compared_res_enrich_list) <-
      paste0("other_", seq_len(length(compared_res_enrich_list)))
  }

  if (!is(compared_res_enrich_list, "list")) {
    stop("You need to provide a list for comparison (even versus one scenario)")
  }

  colnames_res_enrich <- c("gs_id",
                           "gs_description",
                           "gs_pvalue",
                           "gs_genes",
                           "gs_de_count",
                           "gs_bg_count")
  for (i in seq_len(length(compared_res_enrich_list))) {
    this_re <- compared_res_enrich_list[[i]]

    if (!all(colnames_res_enrich %in% colnames(this_re)))
      stop("One of the provided `res_enrich` objects does not respect the format ",
           "required to use in GeneTonic\n",
           "e.g. all required column names have to be present.\n",
           "You might want to use one of the `shake_*` functions to convert it.\n",
           "Required columns: ", paste(colnames_res_enrich, collapse = ", "),
           "\nThis occurred at the element ", i, " in your `compared_res_enrich_list`")

    if (!p_value_column %in% colnames(this_re))
      stop("Required column (p-value) `", p_value_column, "` not found in a component of ",
           "`compared_res_enrich_list` object.",
           "\nThis occurred at the element ", i, " in your `compared_res_enrich_list`")
    if (!color_by %in% colnames(this_re))
      stop("Required column (for coloring) `", color_by, "` not found in a component of ",
           "`compared_res_enrich_list` object.",
           "\nThis occurred at the element ", i, " in your `compared_res_enrich_list`")
  }

  sort_by <- match.arg(sort_by, c("clustered", "first_set"))

  # compared_res_enrich_list
  # append original ref
  all_res_enrichs <- compared_res_enrich_list
  all_res_enrichs[[ref_name]] <- res_enrich

  all_gsids <- lapply(all_res_enrichs, function(arg) arg[["gs_id"]])

  gs_common <- Reduce(intersect, all_gsids)

  if (length(gs_common) == 0) {
    stop("No gene sets have been found in common to the two enrichment results")
  }

  # restrict to the top common n_gs
  gs_common <- gs_common[seq_len(min(n_gs, length(gs_common)))]

  # append scenario info
  res_enrich[["scenario"]] <- ref_name
  compared_res_enrich_list <- lapply(seq_len(length(compared_res_enrich_list)),
                                     function(arg) {
                                       re <- compared_res_enrich_list[[arg]]
                                       re[["scenario"]] <- names(compared_res_enrich_list)[arg]
                                       return(re)
                                     })

  # reduce to common sets
  re_ref <- res_enrich[gs_common, ]
  re_comp <- lapply(seq_len(length(compared_res_enrich_list)),
                    function(arg) {
                      re <- compared_res_enrich_list[[arg]]
                      re <- re[gs_common, ]
                      return(re)
                    })

  merged_res_enh <- rbind(
    re_ref,
    do.call(rbind, re_comp)
  )
  merged_res_enh$logp10 <- -log10(merged_res_enh$gs_pvalue)

  if (sort_by == "first_set") {
    # sorted by category in scenario1
    p <- merged_res_enh %>%
      mutate(gs_description = factor(.data$gs_description, rev(unique(.data$gs_description)))) %>%
      arrange((.data$logp10)) %>%
      ggplot(aes_string(x = "gs_description", y = "logp10")) +
      geom_line(aes_string(group = "scenario", col = "scenario"), size = 3, alpha = 0.7) +
      geom_point(aes_string(fill = "z_score"), size = 4, pch = 21) +
      scale_color_brewer(palette = "Set2") +
      scale_fill_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026") +
      ylim(c(0, NA)) +
      coord_flip() +
      theme_minimal()
  } else if (sort_by == "clustered") {
    # with a nicer sorting - "grouped" by scenario
    nicerorder_terms <- merged_res_enh %>%
      group_by(.data$gs_description) %>%
      mutate(main_category = .data$scenario[which.max(.data$logp10)],
             max_value = max(.data$logp10)) %>%
      arrange(.data$main_category, desc(.data$max_value)) %>%
      dplyr::pull(.data$gs_description)

    p <- merged_res_enh %>%
      # mutate(gs_description=factor(gs_description, unique(gs_description))) %>%
      mutate(gs_description = factor(.data$gs_description, rev(unique(nicerorder_terms)))) %>%
      arrange(desc(.data$logp10)) %>%
      ggplot(aes_string(x = "gs_description", y = "logp10")) +
      geom_line(aes_string(group = "scenario", col = "scenario"), size = 3, alpha = 0.7) +
      scale_color_brewer(palette = "Set2") +
      geom_point(aes_string(fill = "z_score"), size = 4, pch = 21) +
      scale_fill_gradient2(low = "#313695", mid = "#FFFFE5", high = "#A50026") +
      ylim(c(0, NA)) +
      coord_flip() +
      theme_minimal()
  }

  p <- p + labs(x = "Gene set description",
                y = "-log10 p-value",
                col = color_by)

  return(p)
}


#' Plots a heatmap for genes and genesets
#'
#' Plots a heatmap for genes and genesets, useful to spot out intersections across
#' genesets and an overview of them
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
#'
#' @return A `ggplot` object
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
#' gs_summary_heat(res_enrich = res_enrich,
#'                 res_de = res_de,
#'                 annotation_obj = anno_df,
#'                 n_gs = 20)
gs_summary_heat <- function(res_enrich,
                            res_de,
                            annotation_obj,
                            gtl = NULL,
                            n_gs = 80) {
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }

  res_enrich2 <- res_enrich[seq_len(n_gs), ]

  # enriched_gsids <- res_enrich2[["gs_id"]]
  # enriched_gsnames <- res_enrich2[["gs_description"]]
  # enriched_gsdescs <- vapply(enriched_gsids, function(arg) Definition(GOTERM[[arg]]), character(1))

  # rownames(res_enrich) <- enriched_gsids

  gs_expanded <- tidyr::separate_rows(res_enrich2, "gs_genes", sep = ",")
  gs_expanded$log2FoldChange <-
    res_de[annotation_obj$gene_id[match(gs_expanded$gs_genes, annotation_obj$gene_name)], ]$log2FoldChange

  # keep them as factor to prevent rearrangement!
  gs_expanded[["gs_id"]] <- factor(gs_expanded[["gs_id"]], levels = res_enrich2[["gs_id"]])
  gs_expanded[["gs_description"]] <- factor(gs_expanded[["gs_description"]], levels = res_enrich2[["gs_description"]])
  gs_expanded[["gs_genes"]] <- factor(gs_expanded[["gs_genes"]], levels = unique(gs_expanded[["gs_genes"]]))

  p <- ggplot(gs_expanded,
              aes_string(x = "gs_genes", y = "gs_description")) +
    geom_tile(aes_string(fill = "log2FoldChange"),
              color = "white") +
    scale_fill_gradient2(low = muted("deepskyblue"),
                         mid = "lightyellow",
                         high = muted("firebrick"),
                         name = "log2 Fold Change") +
    xlab(NULL) + ylab(NULL) + theme_minimal() +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 75, hjust = 1))

  return(p)
}
