#' Information on a GeneOntology identifier
#'
#' Assembles information, in HTML format, regarding a Gene Ontology identifier
#'
#' Also creates a link to the AmiGO database
#'
#' @param go_id Character, specifying the GeneOntology identifier for which
#' to retrieve information
#' @param res_enrich A `data.frame` object, storing the result of the functional
#' enrichment analysis. If not provided, the experiment-related information is not
#' shown, and only some generic info on the identifier is displayed
#'
#' @return HTML content related to a GeneOntology identifier, to be displayed in
#' web applications (or inserted in Rmd documents)
#' @export
#'
#' @examples
#' go_2_html("GO:0002250")
#' go_2_html("GO:0043368")
go_2_html <- function(go_id, res_enrich = NULL) {
  fullinfo <- GOTERM[[go_id]]
  if (is.null(fullinfo)) {
    return(HTML("GeneOntology term not found!"))
  }

  mycontent <- paste0(
    "<b>GO ID: </b>", .link2amigo(GOID(fullinfo)), "<br>",
    "<b>Term: </b>", Term(fullinfo), "<br></b>",
    ifelse(!is.null(res_enrich),
           paste0("<b>p-value: </b>", res_enrich[(res_enrich$gs_id == go_id), "gs_pvalue"], "</br>",
                  "<b>Z-score: </b>", format(round(res_enrich[(res_enrich$gs_id == go_id), "z_score"], 2), nsmall = 2), "</br>",
                  "<b>Aggregated score: </b>", format(round(res_enrich[(res_enrich$gs_id == go_id), "aggr_score"], 2), nsmall = 2), "</br>",
                  collapse = ""),
           ""),
    "<b>Ontology: </b>", Ontology(fullinfo), "<br><br>",
    "<b>Definition: </b>", Definition(fullinfo), "<br>",
    ## TODO: extra info from the res_enrich?
    paste0(
      unlist(
        lapply(Synonym(fullinfo), function(arg) {
          paste0("<b>Synonym: </b>", arg, "<br>")
        })
      ), collapse = ""
    ),
    ifelse(length(Secondary(fullinfo)) > 0,
           paste0("<b>Secondary: </b>", Secondary(fullinfo), collapse = ""),
           "")
  )
  # ALT IDEA TODO: return as a table?!
  return(HTML(mycontent))

}

#' Link to the AmiGO database
#'
#' @param val A string, with the GO identifier
#'
#' @return HTML for an action button
#' @noRd
.link2amigo <- function(val) {
  sprintf('<a href = "http://amigo.geneontology.org/amigo/term/%s" target = "_blank" class = "btn btn-primary" style = "%s">%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
}


#' Information on a gene
#'
#' Assembles information, in HTML format, regarding a gene symbol identifier
#'
#' Creates links to the NCBI and the GeneCards databases
#'
#' @param gene_id Character specifying the gene identifier for which to retrieve
#' information
#'
#' @return HTML content related to a gene identifier, to be displayed in
#' web applications (or inserted in Rmd documents)
#' @export
#'
#' @examples
#' # TODO
geneinfo_2_html <- function(gene_id) {
  # entrez info
  # genecards?
  # using rentrez?

  mycontent <- paste0(
    "<b>", gene_id, "</b><br>",
    "NCBI link: ",
    .link2ncbi(gene_id), "<br>",
    "GeneCards link: ",
    .link2genecards(gene_id)
  )
  return(HTML(mycontent))
}

#' Link to NCBI database
#'
#' @param val Character, the gene symbol
#'
#' @return HTML for an action button
#' @noRd
.link2ncbi <- function(val) {
  sprintf('<a href = "http://www.ncbi.nlm.nih.gov/gene/?term=%s[sym]" target = "_blank" class = "btn btn-primary" style = "%s">%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
}

#' Link to the GeneCards database
#'
#' @param val Character, the gene symbol of interest
#'
#' @return HTML for an action button
#' @noRd
.link2genecards <- function(val) {
  sprintf('<a href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target = "_blank" class = "btn btn-primary" style = "%s">%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
}

#' Calculate overlap coefficient
#'
#' Calculate similarity coefficient between two sets, based on the overlap
#'
#' @param x Character vector, corresponding to set 1
#' @param y Character vector, set 2
#'
#' @return A numeric value between 0 and 1
#' @export
#'
#' @seealso https://en.wikipedia.org/wiki/Overlap_coefficient
#'
#' @examples
#' # TODO
overlap_coefficient <- function(x, y) {
  length(intersect(x, y)) / min(length(x), length(y))
}

#' Calculate Jaccard Index between two sets
#'
#' Calculate similarity coefficient with the Jaccard Index
#'
#' @param x Character vector, corresponding to set 1
#' @param y Character vector, corresponding to set 2
#'
#' @return A numeric value between 0 and 1
#' @export
#'
#' @examples
#' # TODO
overlap_jaccard_index <- function(x, y) {
  length(intersect(x, y)) / length(unique(c(x, y)))
  # about 2x faster than using union()
}



#' Maps numeric values to color values
#'
#' Maps numeric continuous values to values in a color palette
#'
#' @param x A character vector of numeric values (e.g. log2FoldChange values) to
#' be converted to a vector of colors
#' @param pal A vector of characters specifying the definition of colors for the
#' palette, e.g. obtained via [RColorBrewer::brewer.pal()]
#' @param limits A vector containing the limits of the values to be mapped. If
#' not specified, defaults to the range of values in the `x` vector.
#'
#' @return A vector of colors, each corresponding to an element in the original
#' vector
#' @export
#'
#' @examples
#' # TODO
map2color <- function(x, pal, limits = NULL) {
  if (is.null(limits))
    limits <- range(x)
  pal[findInterval(x, seq(limits[1],
                          limits[2],
                          length.out = length(pal) + 1),
                   all.inside = TRUE)]
}


#' Generate a table from the `DESeq2` results
#'
#' Generate a tidy table with the results of `DESeq2`
#'
#' @param res_de A `DESeqResults` object.
#' @param FDR Numeric value, specifying the significance level for thresholding
#' adjusted p-values. Defaults to NULL, which would return the full set of results
#' without performing any subsetting based on FDR.
#'
#' @return A tidy `data.frame` with the results from differential expression,
#' sorted by adjusted p-value. If FDR is specified, the table contains only genes
#' with adjusted p-value smaller than the value.
#'
#' @export
#'
#' @examples
#' # TODO
deseqresult2df <- function(res_de, FDR = NULL) {
  if (!is(res_de, "DESeqResults"))
    stop("Not a DESeqResults object.")
  res <- as.data.frame(res_de)
  res <- cbind(rownames(res), res)
  names(res)[1] <- "id"
  res$id <- as.character(res$id)
  res <- res[order(res$padj), ]
  if (!is.null(FDR))
    res <- res[!(is.na(res$padj)) & res$padj <= FDR, ]
  res
}

footer <- function() {
  tags$div(
    class = "panel-footer",
    style = "text-align:center",
    tags$div(
      class = "foot-inner",
      list(
        # hr(),
        "GeneTonic is a project developed by Federico Marini in the Bioinformatics division of the ",
        tags$a(href = "http://www.unimedizin-mainz.de/imbei", "IMBEI"),
        "- Institute for Medical Biostatistics, Epidemiology and Informatics", br(),
        "License: ", tags$a(href = "https://opensource.org/licenses/MIT", "MIT"), br(),

        "Development of the GeneTonic package is on ",
        tags$a(href = "https://github.com/federicomarini/GeneTonic", "GitHub")
      )
    )
  )
}


# Shiny resource paths ----------------------------------------------------

.onLoad <- function(libname, pkgname) {
  # Create link to logo
  shiny::addResourcePath("GeneTonic", system.file("www", package = "GeneTonic"))
}


# TODOTODO
# clusterProfiler to common expected format:


# Some constant values ----------------------------------------------------

.actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC; border-color: #2e6da4"
.helpbutton_biocstyle <- "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"

.biocgreen <- "#0092AC"
