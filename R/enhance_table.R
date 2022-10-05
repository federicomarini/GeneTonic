#' Visually enhances a functional enrichment result table
#'
#' Creates a visual summary for the results of a functional enrichment analysis,
#' by displaying also the components of each gene set and their expression change
#' in the contrast of interest
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param res_de  A `DESeqResults` object.
#' @param annotation_obj A `data.frame` object with the feature annotation.
#' information, with at least two columns, `gene_id` and `gene_name`.
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be displayed.
#' @param gs_ids Character vector, containing a subset of `gs_id` as they are
#' available in `res_enrich`. Lists the gene sets to be displayed.
#' @param chars_limit Integer, number of characters to be displayed for each
#' geneset name.
#' @param plot_title Character string, used as title for the plot. If left `NULL`,
#' it defaults to a general description of the plot and of the DE contrast
#'
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
#' enhance_table(res_enrich,
#'   res_de,
#'   anno_df,
#'   n_gs = 10
#' )
enhance_table <- function(res_enrich,
                          res_de,
                          annotation_obj,
                          gtl = NULL,
                          n_gs = 50,
                          gs_ids = NULL,
                          chars_limit = 70,
                          plot_title = NULL) {
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }

  n_gs <- min(n_gs, nrow(res_enrich))

  gs_to_use <- unique(
    c(
      res_enrich$gs_id[seq_len(n_gs)], # the ones from the top
      gs_ids[gs_ids %in% res_enrich$gs_id] # the ones specified from the custom list
    )
  )

  gs_fulllist <- lapply(gs_to_use, function(gs) {
    genes_thisset <- res_enrich[gs, "gs_genes"]
    genes_thisset <- unlist(strsplit(genes_thisset, ","))

    genesid_thisset <- annotation_obj$gene_id[match(genes_thisset, annotation_obj$gene_name)]
    
    # removing the genes not finding a match in the annotation
    no_anno_match <- is.na(genesid_thisset)
    genes_thisset_anno <- genes_thisset[!no_anno_match]
    genesid_thisset_anno <- genesid_thisset[!no_anno_match]
    # ... and informing on which genes might be troublesome
    if (any(no_anno_match)) {
      message("Could not find a match in the annotation for some genes. ",
              "Please inspect your results in detail for geneset ",
              gs,
              " the gene(s) named: ",
              paste0(genes_thisset[no_anno_match], collapse = ", "))
    }
    
    res_thissubset <- res_de[genesid_thisset_anno, ]

    res_thissubset <- as.data.frame(res_thissubset)

    res_thissubset$gene_name <- genes_thisset_anno
    res_thissubset$gs_desc <- as.factor(res_enrich[gs, "gs_description"])
    res_thissubset$gs_id <- res_enrich[gs, "gs_id"]
    # return(as.data.frame(res_thissubset))
    return(res_thissubset)
  })
  gs_fulllist <- do.call(rbind, gs_fulllist)
  # message(dim(gs_fulllist)[1])

  this_contrast <- (sub(".*p-value: (.*)", "\\1", mcols(res_de, use.names = TRUE)["pvalue", "description"]))

  # to have first rows viewed on top
  gs_fulllist <- gs_fulllist[rev(seq_len(nrow(gs_fulllist))), ]
  gs_fulllist$gs_desc <- factor(gs_fulllist$gs_desc, levels = rev(levels(gs_fulllist$gs_desc)))
  max_lfc <- max(abs(range(gs_fulllist$log2FoldChange)))

  p <- ggplot(
    gs_fulllist, aes_string(
      x = "log2FoldChange",
      y = "gs_desc",
      fill = "gs_id",
      text = "gene_name"
    )
  ) +
    scale_x_continuous(limits = c(-max_lfc, max_lfc)) +
    geom_point(alpha = 0.7, shape = 21, size = 2) +
    theme_minimal() +
    geom_vline(aes(xintercept = 0), col = "steelblue", alpha = 0.4) +
    theme(legend.position = "none") +
    scale_y_discrete(
      name = "",
      labels = paste0(
        substr(as.character(unique(gs_fulllist$gs_desc)), 1, chars_limit),
        " | ", unique(gs_fulllist$gs_id)
      )
    ) +
    labs(x = "log2 Fold Change")

  if (is.null(plot_title)) {
    p <- p + ggtitle(paste0("Enrichment overview - ", this_contrast))
  } else {
    p <- p + ggtitle(plot_title)
  }

  return(p)
}


#' Compute aggregated scores for gene sets
#'
#' Computes for each gene set in the `res_enrich` object a Z score and an aggregated
#' score (using the log2FoldChange values, provided in the `res_de`)
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
#' @param aggrfun Specifies the function to use for aggregating the scores for
#' each term. Common values could be `mean` or `median`.
#'
#' @return A `data.frame` with the same columns as provided in the input, with
#' additional information on the `z_score` and the `aggr_score` for each gene set.
#' This information is used by other functions such as [gs_volcano()] or
#' [enrichment_map()]
#'
#' @seealso [gs_volcano()] and [enrichment_map()] make efficient use of the computed
#' aggregated scores
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
#'
#' res_enrich <- get_aggrscores(
#'   res_enrich,
#'   res_de,
#'   anno_df
#' )
get_aggrscores <- function(res_enrich,
                           res_de,
                           annotation_obj,
                           gtl = NULL,
                           aggrfun = mean) {
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }

  gs_expanded <- tidyr::separate_rows(res_enrich, "gs_genes", sep = ",")
  gs_expanded$log2FoldChange <-
    res_de[annotation_obj$gene_id[match(gs_expanded$gs_genes, annotation_obj$gene_name)], ]$log2FoldChange

  gs_aggregated <- lapply(seq_len(nrow(res_enrich)), function(i) {
    this_gsid <- res_enrich$gs_id[i]
    this_subset <- gs_expanded[gs_expanded$gs_id == this_gsid, ]

    upgenes <- sum(this_subset$log2FoldChange > 0)
    downgenes <- sum(this_subset$log2FoldChange < 0)
    z_score <- (upgenes - downgenes) / sqrt(upgenes + downgenes)

    aggr_score <- aggrfun(this_subset$log2FoldChange)
    return(c(
      "DE_count" = nrow(this_subset),
      "Z_score" = z_score,
      "aggr_score" = aggr_score
    ))
  })

  names(gs_aggregated) <- res_enrich$gs_id

  res_enrich$DE_count <- vapply(gs_aggregated, "[", 1, FUN.VALUE = numeric(1))
  res_enrich$z_score <- vapply(gs_aggregated, "[", 2, FUN.VALUE = numeric(1))
  res_enrich$aggr_score <- vapply(gs_aggregated, "[", 3, FUN.VALUE = numeric(1))

  return(res_enrich)
}



#' Distill enrichment results
#'
#' Distill the main topics from the enrichment results, based on the graph derived
#' from constructing an enrichment map
#'
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis.
#' @param res_de A `DESeqResults` object. As for the `dds` parameter, this is
#' also commonly used in the `DESeq2` framework.
#' @param annotation_obj A `data.frame` object, containing two columns, `gene_id`
#' with a set of unambiguous identifiers (e.g. ENSEMBL ids) and `gene_name`,
#' containing e.g. HGNC-based gene symbols.
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#' @param n_gs Integer value, corresponding to the maximal number of gene sets to
#' be used.
#' @param cluster_fun Character, referring to the name of the function used for
#' the community detection in the enrichment map graph. Could be one of "cluster_markov",
#' "cluster_louvain", or "cluster_walktrap", as they all return a `communities`
#' object.
#'
#' @return A list containing three objects:
#' - the distilled table of enrichment, `distilled_table`, where the new meta-genesets
#' are identified and defined, specifying e.g. the names of each component, and the
#' genes associated to these.
#' - the distilled graph for the enrichment map, `distilled_em`, with the information
#' on the membership
#' - the original `res_enrich`, augmented with the information of the membership
#' related to the meta-genesets
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
#' distilled <- distill_enrichment(res_enrich,
#'   res_de,
#'   annotation_obj,
#'   n_gs = 100,
#'   cluster_fun = "cluster_markov"
#' )
#' colnames(distilled$distilled_table)
#' distilled$distilled_em
distill_enrichment <- function(res_enrich,
                               res_de,
                               annotation_obj,
                               gtl = NULL,
                               n_gs = nrow(res_enrich),
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

  n_gs <- min(n_gs, nrow(res_enrich))

  em <- enrichment_map(res_enrich,
    res_de,
    annotation_obj,
    n_gs = n_gs
  )

  # subset accordingly
  res_enrich <- res_enrich[seq_len(n_gs), ]

  gs_communities <- cluster_fun(em)
  res_enrich$gs_membership <- factor(gs_communities$membership)
  V(em)$membership <- gs_communities$membership
  V(em)$color <- gs_communities$membership


  # aggregate the results according to the defined gs_membership column

  distilled_res <- data.frame(
    metags_cluster = unique(res_enrich$gs_membership),
    metags_n_gs = NA,
    metags_genes = NA,
    metags_n_genes = NA,
    metags_gsidlist = NA,
    metags_gsdesclist = NA,
    metags_msgs = NA,
    metags_mcgs = NA,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(distilled_res$metags_cluster)) {
    # message(i)
    # message(distilled_res$metags_cluster[i])
    subset_enrich <- res_enrich[res_enrich$gs_membership == distilled_res$metags_cluster[i], ]
    distilled_res[i, "metags_n_gs"] <- nrow(subset_enrich)

    all_genes_singlevec <- unique(strsplit(paste0(subset_enrich$gs_genes, collapse = ","), ",")[[1]])
    distilled_res[i, "metags_genes"] <- paste0(all_genes_singlevec, collapse = ",")
    distilled_res[i, "metags_n_genes"] <- length(all_genes_singlevec)

    distilled_res[i, "metags_gsidlist"] <- paste0(subset_enrich$gs_id, collapse = ",") # nested list or collapsed?
    distilled_res[i, "metags_gsdesclist"] <- paste0(subset_enrich$gs_description, collapse = ",")

    most_sig <- which.min(subset_enrich$gs_pvalue)
    distilled_res[i, "metags_msgs"] <- paste0(
      subset_enrich$gs_id[most_sig], "|",
      subset_enrich$gs_description[most_sig], "|",
      subset_enrich$gs_pvalue[most_sig]
    )

    mgs_graph <- induced_subgraph(em, subset_enrich$gs_description)
    distilled_res[i, "metags_mcgs"] <- names(which.max(strength(mgs_graph)))
  }
  # later add something maybe even based on NLP/wordcloud or so


  # lapply(unique(res_enrich$gs_membership), function(gs_cluster) {
  # cur_clus <- gs_cluster
  # message(cur_clus)
  # clu_n_gs <- res_enrich$gs_membership
  # })

  return(
    list(
      distilled_table = distilled_res,
      distilled_em = em,
      res_enrich = res_enrich
    )
  )
}
