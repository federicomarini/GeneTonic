#' Title TODO
#'
#' TODO
#'
#' @param go_id _TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
go_2_html <- function(go_id) {
  fullinfo <- GOTERM[[go_id]]
  if(is.null(fullinfo)) {
    return(HTML("GeneOntology term not found!"))
  }

  mycontent <- paste0(
    "<b>GO ID: </b>", .link2amigo(GOID(fullinfo)), "<br>",
    "<b>Term: </b>", Term(fullinfo),"<br></b>",
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

.link2amigo <- function(val) {
  sprintf('<a href = "http://amigo.geneontology.org/amigo/term/%s" target = "_blank" class = "btn btn-primary" style = "%s">%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
}


#' Title TODO
#'
#' TODO
#'
#' @param gene_id TODO
#'
#' @return TODO
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

.link2ncbi <- function(val) {
  sprintf('<a href = "http://www.ncbi.nlm.nih.gov/gene/?term=%s[sym]" target = "_blank" class = "btn btn-primary" style = "%s">%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
}

.link2genecards <- function(val) {
  sprintf('<a href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target = "_blank" class = "btn btn-primary" style = "%s">%s</a>',
          val,
          .actionbutton_biocstyle,
          val)
}

#' Title
#'
#' TODO
#'
#' @param x TODO
#' @param y TODO
#'
#' @return TODO
#' @export
#'
#' @seealso https://en.wikipedia.org/wiki/Overlap_coefficient
#'
#' @examples
#' # TODO
overlap_coefficient <- function(x, y) {
  length(intersect(x, y)) / min(length(x), length(y))
}

#' Title
#'
#' TODO
#'
#' @param x TODO
#' @param y TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
overlap_jaccard_index <- function(x, y) {
  length(intersect(x, y)) / length(unique(c(x, y)))
  # about 2x faster than using union()
}



#' Title
#'
#' TODO
#'
#' @param x TODO
#' @param pal TODO
#' @param limits TODO
#'
#' @return TODO
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


#' Title TODO
#'
#' TODO
#'
#' @param deseqresult TODO
#' @param FDR TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
deseqresult2df <- function(deseqresult, FDR = NULL) {
  if (!is(deseqresult, "DESeqResults"))
    stop("Not a DESeqResults object.")
  res <- as.data.frame(deseqresult)
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

# TODOTODO
# clusterProfiler to common expected format:


# Some constant values ----------------------------------------------------

.actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC; border-color: #2e6da4"

.biocgreen <- "#0092AC"


