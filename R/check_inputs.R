#' Checking the input objects for GeneTonic
#'
#' Checking the input objects for GeneTonic, whether these are all set for running
#' the app
#'
#' Some suggestions on the requirements for each parameter are returned in the
#' error messages.
#'
#' @param dds A `DESeqDataSet` object, normally obtained after running your data
#' through the `DESeq2` framework.
#' @param res_de A `DESeqResults` object. As for the `dds` parameter, this is
#' also commonly used in the `DESeq2` framework.
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
#' @param annotation_obj A `data.frame` object, containing two columns, `gene_id`
#' with a set of unambiguous identifiers (e.g. ENSEMBL ids) and `gene_name`,
#' containing e.g. HGNC-based gene symbols.
#' @param verbose Logical, to control level of verbosity of the messages generated
#' @return Invisible NULL
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
#' checkup_GeneTonic(dds = dds_macrophage,
#'                   res_de = res_de,
#'                   res_enrich = res_enrich,
#'                   annotation_obj = anno_df)
#' # if all is fine, it should return an invisible NULL and a simple message
checkup_GeneTonic <- function(dds,
                              res_de,
                              res_enrich,
                              annotation_obj,
                              verbose = FALSE) {

  # checking dds
  if (!is(dds, "DESeqDataSet"))
    stop("The provided `dds` is not a DESeqDataSet object, please check your input parameters")

  # check normalization factors/size factors are in
  if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
    warning("Could not find size factors for the provided `dds` object...",
            "\nPlease compute them first by calling\n",
            "'dds <- estimateSizeFactors(dds)'")
    # dds <- estimateSizeFactors(dds)
  }

  # checking res_de
  if (!is(res_de, "DESeqResults"))
    stop("The provided `res_de` is not a DESeqResults object, please check your input parameters")

  # checking that results and dds are related
  ## at least a subset of dds should be in res
  if (!all(rownames(res_de) %in% rownames(dds)))
    warning("It is likely that the provided `dds` and `res_de` objects are not related ",
            "to the same dataset (the row names of the results are not all in the dds). ",
            "Are you sure you want to proceed?")

  # checking res_enrich
  colnames_res_enrich <- c("gs_id",
                           "gs_description",
                           "gs_pvalue",
                           "gs_genes",
                           "gs_de_count",
                           "gs_bg_count") # gs_ontology? more?
  if (!all(colnames_res_enrich %in% colnames(res_enrich)))
    stop("The provided `res_enrich` object does not respect the format required to use in GeneTonic\n",
         "e.g. all required column names have to be present.\n",
         "You might want to use one of the `shaker_*` functions to convert it.\n",
         "Required columns: ", paste(colnames_res_enrich, collapse = ", "))

  # checking annotation_obj
  colnames_annotation_obj <- c("gene_id",
                               "gene_name")
  if (!all(colnames_annotation_obj %in% colnames(annotation_obj)))
    stop("The provided `annotation_obj` object does not respect the format required to use in GeneTonic\n",
         "e.g. all required column names have to be present.\n",
         "You can use e.g. `pcaExplorer::get_annotation_orgdb()` for this purpose.\n",
         "Required columns: ", paste(colnames_annotation_obj, collapse = ", "))

  if (verbose)
    message("All set to enjoy GeneTonic!")
  invisible(NULL)
}


#' Checking the `gtl` input object for GeneTonic
#'
#' Checking the `gtl` ("GeneTonic list") input object for GeneTonic, with the 
#' correct content and format expected
#'
#' Some suggestions on the requirements for the `gtl` are returned in the
#' error messages.
#'
#' @param gtl A `DESeqDataSet` object, normally obtained after running your data
#' through the `DESeq2` framework.
#' This list should contain
#' - in the `dds` slot: A `DESeqDataSet` object
#' - in the `res_de`: A `DESeqResults` object
#' - in the `res_enrich`: A `data.frame` object, storing the result of the 
#' functional enrichment analysis
#' - in the `annotation_obj`: A `data.frame` object, containing two columns, 
#' `gene_id` with a set of unambiguous identifiers (e.g. ENSEMBL ids) and 
#' `gene_name`, containing e.g. HGNC-based gene symbols.
#' @param verbose Logical, to control level of verbosity of the messages generated
#'
#' @return Invisible NULL
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
#' gtl <- list(dds = dds_macrophage,
#'             res_de = res_de,
#'             res_enrich = res_enrich,
#'             annotation_obj = anno_df)
#' 
#' checkup_gtl(gtl)
#' # if all is fine, it should return an invisible NULL and a simple message
checkup_gtl <- function(gtl,
                        verbose = FALSE) {
  
  # checking gtl
  if (!is(gtl, "list"))
    stop("The provided `gtl` is not a list object, please check your input parameters")
  
  if(!all(c("dds", "res_de", "res_enrich", "annotation_obj") %in% names(gtl)))
    stop("The names of the provided gtl object are not specified as expected")
  
  dds <- gtl$dds
  res_de <- gtl$res_de
  res_enrich <- gtl$res_enrich
  annotation_obj <- gtl$annotation_obj
  
  if (verbose)
    message("Parameter passed as a list...")
  checkup_GeneTonic(dds, res_de, res_enrich, annotation_obj, verbose)
  
  invisible(NULL)
}


