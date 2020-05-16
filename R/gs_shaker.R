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
#' @param obj An `enrichResult` object, obtained via `clusterProfiler`
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
  if (!all(colnames(my_david) %in% exp_colnames))
    warning("I could not find some of the usual column names from the DAVID output exported to file")
  
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
