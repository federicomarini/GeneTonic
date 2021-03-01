# for all the available enrichment results type, we need
  # - gs_id
  # - gs_description
  # - gs_pvalue
  # - gs_genes
 # optional:
  # - gs_de_count
  # - gs_bg_count
# - plus: cleverly the row name should be gs_id (consistent with clever indexing)
# - plus: for practical reasons, I could use the original columns - or store it in a slot?


#' Convert an enrichResult object
#'
#' Convert an enrichResult object for straightforward use in [GeneTonic()]
#'
#' This function is able to handle the output of `clusterProfiler` and `reactomePA`,
#' as they both return an object of class `enrichResult` - and this in turn
#' contains the information required to create correctly a `res_enrich` object.
#'
#' @param obj An `enrichResult` object, obtained via `clusterProfiler` (or also
#' via `reactomePA`)
#'
#' @return A `data.frame` compatible for use in [GeneTonic()] as `res_enrich`
#' @export
#'
#' @family shakers
#'
#' @examples
#' # dds
#' library("macrophage")
#' library("DESeq2")
#' data(gse)
#' dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#'
#' # res object
#' data(res_de_macrophage, package = "GeneTonic")
#' res_de <- res_macrophage_IFNg_vs_naive
#' de_symbols_IFNg_vs_naive <- res_macrophage_IFNg_vs_naive[
#'   (!(is.na(res_macrophage_IFNg_vs_naive$padj))) &
#'   (res_macrophage_IFNg_vs_naive$padj <= 0.05), "SYMBOL"]
#' bg_ids <- rowData(dds_macrophage)$SYMBOL[rowSums(counts(dds_macrophage)) > 0]
#' \dontrun{
#' library("clusterProfiler")
#' library("org.Hs.eg.db")
#' ego_IFNg_vs_naive <- enrichGO(gene = de_symbols_IFNg_vs_naive,
#'                               universe      = bg_ids,
#'                               keyType       = "SYMBOL",
#'                               OrgDb         = org.Hs.eg.db,
#'                               ont           = "BP",
#'                               pAdjustMethod = "BH",
#'                               pvalueCutoff  = 0.01,
#'                               qvalueCutoff  = 0.05,
#'                               readable      = FALSE)
#'
#' res_enrich <- shake_enrichResult(ego_IFNg_vs_naive)
#' head(res_enrich)
#' }
shake_enrichResult <- function(obj) {
  if (!is(obj, "enrichResult"))
    stop("Provided object must be of class `enrichResult`")

  if (is.null(obj@result$geneID)) {
    stop("You are providing an object where the gene symbols are not specified, ",
         "this is required for running GeneTonic properly.")
  }

  message("Found ", nrow(obj@result), " gene sets in `enrichResult` object, of which ", nrow(as.data.frame(obj)), " are significant.")
  message("Converting for usage in GeneTonic...")

  fullresults <- obj@result

  mydf <- data.frame(
    gs_id = fullresults$ID,
    gs_description = fullresults$Description,
    gs_pvalue = fullresults$pvalue,
    gs_genes = gsub("/", ",", fullresults$geneID),
    gs_de_count = fullresults$Count,
    gs_bg_count = unlist(lapply(strsplit(fullresults$BgRatio, "/"), function(arg) arg[[1]])),
    gs_ontology = obj@ontology,
    GeneRatio = fullresults$GeneRatio,
    BgRatio = fullresults$BgRatio,
    p.adjust = fullresults$p.adjust,
    qvalue = fullresults$qvalue,
    stringsAsFactors = FALSE
  )

  rownames(mydf) <- mydf$gs_id

  return(mydf)
}




#' Convert a topGOtableResult object
#'
#' Convert a topGOtableResult object for straightforward use in [GeneTonic()]
#'

#'
#' @param obj A `topGOtableResult` object
#' @param p_value_column Character, specifying which column the p value for
#' enrichment has to be used. Example values are "p.value_elim" or "p.value_classic"
#'
#' @return A `data.frame` compatible for use in [GeneTonic()] as `res_enrich`
#' @export
#'
#' @family shakers
#'
#' @examples
#'
#' # res_enrich object
#' data(res_enrich_macrophage, package = "GeneTonic")
#'
#' res_enrich <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)
#'
shake_topGOtableResult <- function(obj,
                                   p_value_column = "p.value_elim") {

  if(!all(c("GO.ID", "Term", "Annotated", "Significant", "Expected", "p.value_classic") %in%
          colnames(obj))) {
    stop("The provided object must be of in the format specified by the `pcaExplorer::topGOtable` function")
  }

  if(!p_value_column %in% colnames(obj)) {
    stop("You specified a column for the p-value which is not contained in the provided object. \n",
         "Please check the colnames of your object in advance.")
  }

  if(!"genes" %in% colnames(obj)) {
    stop("The column `genes` is not present in the provided object and is required for properly running GeneTonic.",
         "\nMaybe you did set `addGeneToTerms` to FALSE in the call to `pcaExplorer::topGOtable`?")
  }

  # Thought: store somewhere the ontology if possible - in an extra column?
  message("Found ", nrow(obj), " gene sets in `topGOtableResult` object.")
  message("Converting for usage in GeneTonic...")

  fullresults <- obj

  mydf <- data.frame(
    gs_id = fullresults$GO.ID,
    gs_description = fullresults$Term,
    gs_pvalue = fullresults[[p_value_column]],
    gs_genes = fullresults$genes,
    gs_de_count = fullresults$Significant,
    gs_bg_count = fullresults$Annotated,
    # gs_ontology = obj@ontology,
    Expected = fullresults$Expected,
    stringsAsFactors = FALSE
  )

  rownames(mydf) <- mydf$gs_id

  return(mydf)
}

## potential extensions, if the formats/classes get defined in a stable/compatible manner
# shake_goseqResult?
# shake_viseago?
# shake_GSECA?


#' Convert the output of DAVID
#'
#' Convert the output of DAVID for straightforward use in [GeneTonic()]
#'
#' @param david_output_file The location of the text file output, as exported from
#' DAVID
#'
#' @return A `data.frame` compatible for use in [GeneTonic()] as `res_enrich`
#' @export
#'
#' @family shakers
#'
#' @examples
#' david_output_file <- system.file("extdata",
#'                                  "david_output_chart_BPonly_ifng_vs_naive.txt",
#'                                  package = "GeneTonic")
#' res_enrich <- shake_davidResult(david_output_file)
shake_davidResult <- function(david_output_file) {

  if(!file.exists(david_output_file))
    stop("File not found")

  my_david <- read.delim(david_output_file, header = TRUE, sep = "\t")
  # careful, names are auto-sanitized

  exp_colnames <- c("Category", "Term", "Count", "X.", "PValue", "Genes",
                    "List.Total", "Pop.Hits", "Pop.Total", "Fold.Enrichment",
                    "Bonferroni", "Benjamini", "FDR")
  if (!all(exp_colnames %in% colnames(my_david) ))
    stop("I could not find some of the usual column names from the DAVID output exported to file")

  message("Found ", nrow(my_david), " gene sets in the file output from DAVID of which ", sum(my_david$PValue <= 0.05), " are significant (p-value <= 0.05).")
  message("Converting for usage in GeneTonic...")

  mydf <- data.frame(
    gs_id = unlist(lapply(strsplit(my_david$Term, "~"), function(arg) arg[[1]])),
    gs_description = unlist(lapply(strsplit(my_david$Term, "~"), function(arg) arg[[2]])),
    gs_pvalue = my_david$PValue,
    gs_genes = gsub(", ", ",", my_david$Genes),
    gs_de_count = my_david$Count,
    gs_bg_count = my_david$Pop.Hits,
    gs_ontology = my_david$Category,
    gs_generatio = my_david$Count/my_david$List.Total,
    gs_bgratio = my_david$Pop.Hits/my_david$Pop.Total,
    gs_foldenrich  = my_david$Fold.Enrichment,
    gs_bonferroni = my_david$Bonferroni,
    gs_benjamini = my_david$Benjamini,
    gs_FDR = my_david$FDR,
    stringsAsFactors = FALSE
  )

  rownames(mydf) <- mydf$gs_id

  return(mydf)
}



#' Convert the output of Enrichr
#'
#' Convert the output of Enrichr for straightforward use in [GeneTonic()]
#'
#' @param enrichr_output_file The location of the text file output, as exported from
#' Enrichr
#' @param enrichr_output A data.frame with the output of `enrichr`, related to a
#' specific set of genesets. Usually it is one of the members of the list returned
#' by the initial call to `enrichr`.
#'
#' @return A `data.frame` compatible for use in [GeneTonic()] as `res_enrich`
#' @export
#'
#' @family shakers
#'
#' @examples
#' # library("enrichR")
#' # dbs <- c("GO_Molecular_Function_2018",
#' #          "GO_Cellular_Component_2018",
#' #          "GO_Biological_Process_2018",
#' #          "KEGG_2019_Human",
#' #          "Reactome_2016",
#' #          "WikiPathways_2019_Human")
#' # degenes <- (deseqresult2df(res_macrophage_IFNg_vs_naive, FDR = 0.01)$SYMBOL)
#' # if called directly within R...
#' # enrichr_output_macrophage <- enrichr(degenes, dbs)
#' # or alternatively, if downloaded from the website in tabular format
#' enrichr_output_file <- system.file("extdata",
#'                                    "enrichr_tblexport_IFNg_vs_naive.txt",
#'                                    package = "GeneTonic")
#' res_from_enrichr <- shake_enrichrResult(enrichr_output_file = enrichr_output_file)
#' # res_from_enrichr2 <- shake_enrichrResult(
#' #   enrichr_output = enrichr_output_macrophage[["GO_Biological_Process_2018"]])
shake_enrichrResult <- function(enrichr_output_file,
                                enrichr_output = NULL) {
  if (is.null(enrichr_output)) {
    if(!file.exists(enrichr_output_file))
      stop("File not found")
    enrichr_output <- read.delim(enrichr_output_file, header = TRUE, sep = "\t")
  }

  if (!is.null(enrichr_output)) {
    # if still a list, might need to select the appropriate element
    if (is(enrichr_output, "list"))
      stop("Expecting a data.frame object. Maybe you are providing the list",
           " containing it? You could do so by selecting the appropriate element",
           " of the list")
  }

  exp_colnames <- c("Term", "Overlap", "P.value", "Adjusted.P.value",
                    "Old.P.value", "Old.Adjusted.P.value", "Odds.Ratio",
                    "Combined.Score", "Genes")
  if (!all(colnames(enrichr_output) %in% exp_colnames))
    stop("I could not find some of the usual column names from the Enrichr output")

  message("Found ", nrow(enrichr_output), " gene sets in the file output from Enrichr of which ", sum(enrichr_output$P.value <= 0.05), " are significant (p-value <= 0.05).")
  message("Converting for usage in GeneTonic...")

  # TODO: split up id and term - or just keep em the same?!
  # this does work for go term as they encode it...
  mydf <- data.frame(
    gs_id = gsub("\\)", "", gsub("^.* \\(", "", enrichr_output$Term)),
    gs_description = gsub(" \\(GO.*$", "", enrichr_output$Term),
    gs_pvalue = enrichr_output$P.value,
    gs_genes = gsub(";", ",", enrichr_output$Genes),
    gs_de_count = as.numeric(
      unlist(lapply(strsplit(enrichr_output$Overlap, "/"), function(arg) arg[[1]]))),
    gs_bg_count = as.numeric(
      unlist(lapply(strsplit(enrichr_output$Overlap, "/"), function(arg) arg[[2]]))),
    gs_adj_pvalue = enrichr_output$Adjusted.P.value,
    stringsAsFactors = FALSE
  )

  rownames(mydf) <- mydf$gs_id

  return(mydf)
}


#' Convert the output of g:Profiler
#'
#' Convert the output of g:Profiler for straightforward use in [GeneTonic()]
#'
#' @param gprofiler_output_file The location of the text file output, as exported from
#' g:Profiler
#' @param gprofiler_output A data.frame with the output of `gost()` in `gprofiler2`.
#' Usually it is one of the members of the list returned by the initial call to `gost`.
#'
#' @return A `data.frame` compatible for use in [GeneTonic()] as `res_enrich`
#' @export
#'
#' @family shakers
#'
#' @examples
#' # degenes <- (deseqresult2df(res_macrophage_IFNg_vs_naive, FDR = 0.01)$SYMBOL)
#' # if called directly withín R...
#' # enrichr_output_macrophage <- enrichr(degenes, dbs)
#' # or alternatively, if downloaded from the website in tabular format
#' gprofiler_output_file <- system.file(
#'   "extdata",
#'   "gProfiler_hsapiens_5-25-2020_tblexport_IFNg_vs_naive.csv",
#'   package = "GeneTonic")
#' res_from_gprofiler <- shake_gprofilerResult(gprofiler_output_file = gprofiler_output_file)
#'
#' data(gostres_macrophage, package = "GeneTonic")
#' res_from_gprofiler_2 <- shake_gprofilerResult(
#'   gprofiler_output = gostres_macrophage$result
#' )
shake_gprofilerResult <- function(gprofiler_output_file,
                                  gprofiler_output = NULL) {
  if (is.null(gprofiler_output)) {
    # expecting text input
    if(!file.exists(gprofiler_output_file))
      stop("File not found")
    gprofiler_output <- read.delim(gprofiler_output_file, header = TRUE, sep = ",")

    exp_colnames_textual <- c("source", "term_name", "term_id", "adjusted_p_value",
                              "negative_log10_of_adjusted_p_value", "term_size",
                              "query_size", "intersection_size", "effective_domain_size",
                              "intersections")

    if (!all(exp_colnames_textual %in% colnames(gprofiler_output)))
      stop("I could not find some of the usual column names from the g:Profiler output.",
           " A possible reason could be that you did not specify `evcodes = TRUE`?",
           " This is required to fill in all the required fields of `res_enrich`")

    message("Found ", nrow(gprofiler_output), " gene sets in the file output from Enrichr of which ", sum(gprofiler_output$adjusted_p_value <= 0.05), " are significant (p-value <= 0.05).")
    message("Converting for usage in GeneTonic...")

    mydf <- data.frame(
      gs_id = gprofiler_output$term_id,
      gs_description = gprofiler_output$term_name,
      gs_pvalue = gprofiler_output$adjusted_p_value,
      gs_genes = gprofiler_output$intersections,
      gs_de_count = gprofiler_output$intersection_size,
      gs_bg_count = gprofiler_output$term_size,
      gs_adj_pvalue = gprofiler_output$adjusted_p_value,
      stringsAsFactors = FALSE
    )

  } else {
    # using directly the output from the call from gprofiler2
    # if still a list, might need to select the appropriate element
    if (is(gprofiler_output, "list"))
      stop("Expecting a data.frame object. Maybe you are providing the list",
           " containing it? You could do so by selecting the appropriate element",
           " of the list")

    exp_colnames_rcall <- c("query", "significant", "p_value", "term_size", "query_size",
                            "intersection_size", "precision", "recall",
                            "term_id", "source", "term_name", "effective_domain_size",
                            "source_order", "parents", "evidence_codes", "intersection")

    if (!all(colnames(gprofiler_output) %in% exp_colnames_rcall))
      stop("I could not find some of the usual column names from the g:Profiler output.",
           " A possible reason could be that you did not specify `evcodes = TRUE`?",
           " This is required to fill in all the required fields of `res_enrich`")

    message("Found ", nrow(gprofiler_output), " gene sets in the file output from Enrichr of which ", sum(gprofiler_output$p_value <= 0.05), " are significant (p-value <= 0.05).")
    message("Converting for usage in GeneTonic...")

    mydf <- data.frame(
      gs_id = gprofiler_output$term_id,
      gs_description = gprofiler_output$term_name,
      gs_pvalue = gprofiler_output$p_value,
      gs_genes = gprofiler_output$intersection,
      gs_de_count = gprofiler_output$intersection_size,
      gs_bg_count = gprofiler_output$term_size,
      gs_adj_pvalue = gprofiler_output$p_value,
      gs_ontology = gprofiler_output$source,
      stringsAsFactors = FALSE
    )
  }

  rownames(mydf) <- mydf$gs_id

  return(mydf)
}



#' Convert the output of fgsea
#'
#' Convert the output of fgsea for straightforward use in [GeneTonic()]
#'
#' @param fgsea_output A data.frame with the output of `fgsea()` in `fgsea`.
#'
#' @return A `data.frame` compatible for use in [GeneTonic()] as `res_enrich`
#' @export
#'
#' @family shakers
#'
#' @examples
#' data(fgseaRes, package = "GeneTonic")
#' res_from_fgsea <- shake_fgseaResult(fgseaRes)
shake_fgseaResult <- function(fgsea_output) {
  if (!is(fgsea_output, "data.frame")) {
    stop("fgsea output should be a data.frame!")
  }
  exp_colnames <- c("pathway", "pval", "padj", "ES", "NES",
                    "nMoreExtreme", "size", "leadingEdge")
  if (!all(colnames(fgsea_output) %in% exp_colnames))
    stop("I could not find some of the usual column names from the fgsea output.",
         " Maybe you performed additional processing/filtering steps?")

  if (!is(fgsea_output$leadingEdge, "list")) {
    stop("Expecting 'leadingEdge' column to be a list")
  }

  message("Found ", nrow(fgsea_output), " gene sets in the file output from Enrichr of which ", sum(fgsea_output$padj <= 0.05), " are significant (p-value <= 0.05).")
  message("Converting for usage in GeneTonic...")

  message(
    "Using the content of the 'leadingEdge' column to generate the 'gs_genes' for GeneTonic...",
    " If you have that information available directly, please adjust the content accordingly.",
    "\n\nUsing the set of the leadingEdge size to compute the 'gs_de_count'"
  )

  message("\n\nfgsea is commonly returning no identifier for the gene sets used.",
          " Please consider editing the 'gs_id' field manually according to the gene set you",
          " provided")

  mydf <- data.frame(
    gs_id = fgsea_output$pathway,
    gs_description = fgsea_output$pathway,
    gs_pvalue = fgsea_output$pval,
    gs_genes = vapply(fgsea_output$leadingEdge,
                      function (arg) paste(arg, collapse = ","), character(1)),
    gs_de_count = lengths(fgsea_output$leadingEdge),
    gs_bg_count = fgsea_output$size,
    gs_NES = fgsea_output$NES,
    gs_adj_pvalue = fgsea_output$padj,
    stringsAsFactors = FALSE
  )

  rownames(mydf) <- mydf$gs_id

  # consider re-sorting by p-value?


  return(mydf)
}


