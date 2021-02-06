#' export_for_iSEE
#' 
#' Combine data from a typical DESeq2 run
#' 
#' Combines the DESeqDataSet input and DESeqResults into a SummarizedExperiment
#' object, which can be readily explored with iSEE.
#' 
#' A typical usage would be after running the DESeq2 pipeline and/or after exploring
#' the functional enrichment results with [GeneTonic()]
#'
#' @param dds A \code{\link{DESeqDataSet}} object.
#' @param res_de A \code{\link{DESeqResults}} object.
#' @param gtl A `GeneTonic`-list object, containing in its slots the arguments
#' specified above: `dds`, `res_de`, `res_enrich`, and `annotation_obj` - the names
#' of the list _must_ be specified following the content they are expecting
#' 
#' @return A `SummarizedExperiment` object, with raw counts, normalized counts, and 
#' variance-stabilizing transformed counts in the assay slots; and with colData 
#' and rowData extracted from the corresponding input parameters - mainly the 
#' results for differential expression analysis.
#' 
#' @export
#'
#' @examples
#' library("macrophage")
#' library("DESeq2")
#'
#' # dds object
#' data("gse", package = "macrophage")
#' dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' dds_macrophage <- estimateSizeFactors(dds_macrophage)
#'
#' # res object
#' data(res_de_macrophage, package = "GeneTonic")
#' res_de <- res_macrophage_IFNg_vs_naive
#'
#' # now everything is in place to launch the app
#' # dds_macrophage <- DESeq2::DESeq(dds_macrophage)
#' se_macrophage <- export_for_iSEE(dds_macrophage, res_de)
#' # iSEE(se_macrophage)
export_for_iSEE <- function(dds, 
                            res_de,
                            gtl = NULL) {
  
  if (!is.null(gtl)) {
    checkup_gtl(gtl)
    dds <- gtl$dds
    res_de <- gtl$res_de
    res_enrich <- gtl$res_enrich
    annotation_obj <- gtl$annotation_obj
  }
  
  # sanity checks on the objects
  if (length(setdiff(rownames(dds), rownames(res_de))) != 0 | 
      length(setdiff(rownames(res_de), rownames(dds))) != 0) {
    message("Attempting to combine a dds and a res_de object where not all",
            "identifiers are shared. Subsetting to the intersection of these..." )
  }
  
  keep <- intersect(rownames(dds), rownames(res_de))
  dds <- dds[keep, ]
  res_de <- res_de[keep, ]
  
  # dds to vst
  vst <- vst(dds)
  
  # initialize the container
  se <- SummarizedExperiment(
    assays = list(
      counts = counts(dds),
      normcounts = counts(dds,normalized = TRUE),
      vst_counts = assay(vst)
    )
  )
  
  # adding colData, taken directly from the DESeqDataSet object
  colData(se) <- colData(dds)
  
  # extract contrast info
  this_contrast <- sub(".*p-value: (.*)","\\1",mcols(res_de, use.names=TRUE)["pvalue","description"])
  
  # getting the rowData from the dds itself
  rdd <- rowData(dds)
  
  # modifying in advance the DESeqResults object
  res_de$log10_baseMean <- log10(res_de$baseMean)
  res_de$log10_pvalue <- -log10(res_de$pvalue)
  
  if ("dispersion" %in% colnames(rdd)) {
    # and for the rowData
    rdd$log10_dispersion <- log10(rdd$dispersion)
  }  
  
  # adding rowData to se
  rowData(se)[[paste0("DESeq2_",gsub(" ","_",this_contrast))]] <- res_de
  
  # merging in the existing rowData slot
  rowData(se) <- cbind(rowData(se), rdd)
  
  return(se)
}
