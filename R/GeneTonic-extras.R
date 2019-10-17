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


.actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC; border-color: #2e6da4"

.biocgreen <- "#0092AC"


