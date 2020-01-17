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
#' shown, and only some generic info on the identifier is displayed.
#' See more in the main function, [GeneTonic()], to check the
#' formatting requirements (a minimal set of columns should be present).
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
#' geneinfo_2_html("ACTB")
#' geneinfo_2_html("Pf4")
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
#' a <- seq(1, 21, 2)
#' b <- seq(1, 11, 2)
#' overlap_coefficient(a,b)
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
#' a <- seq(1, 21, 2)
#' b <- seq(1, 11, 2)
#' overlap_jaccard_index(a,b)
overlap_jaccard_index <- function(x, y) {
  length(intersect(x, y)) / length(unique(c(x, y)))
  # about 2x faster than using union()
}






#' Style DT color bars
#'
#' Style DT color bars for values that diverge from 0.
#'
#' @details This function draws background color bars behind table cells in a column,
#' width the width of bars being proportional to the column values *and* the color
#' dependent on the sign of the value.
#'
#' A typical usage is for values such as `log2FoldChange` for tables resulting from
#' differential expression analysis.
#' Still, the functionality of this can be quickly generalized to other cases -
#' see in the examples.
#'
#' The code of this function is heavily inspired from styleColorBar, and borrows
#' at full hands from an excellent post on StackOverflow -
#' https://stackoverflow.com/questions/33521828/stylecolorbar-center-and-shift-left-right-dependent-on-sign/33524422#33524422
#'
#' @param data The numeric vector whose range will be used for scaling the table
#' data from 0-100 before being represented as color bars. A vector of length 2
#' is acceptable here for specifying a range possibly wider or narrower than the
#' range of the table data itself.
#' @param color_pos The color of the bars for the positive values
#' @param color_neg The color of the bars for the negative values
#'
#' @return This function generates JavaScript and CSS code from the values
#' specified in R, to be used in DT tables formatting.
#'
#' @export
#'
#' @examples
#'
#' data(res_de_macrophage, package = "GeneTonic")
#' res_df <- deseqresult2df(res_macrophage_IFNg_vs_naive)
#' library("magrittr")
#' library("DT")
#' DT::datatable(res_df [1:50, ],
#'               options = list(
#'                 pageLength = 25,
#'                 columnDefs = list(
#'                   list(className = "dt-center", targets = "_all")
#'                 )
#'               )
#' ) %>%
#'   formatRound(columns = c("log2FoldChange"), digits = 3) %>%
#'   formatStyle(
#'     "log2FoldChange",
#'     background = styleColorBar_divergent(res_df$log2FoldChange,
#'                                          scales::alpha("navyblue", 0.4),
#'                                          scales::alpha("darkred", 0.4)),
#'     backgroundSize = "100% 90%",
#'     backgroundRepeat = "no-repeat",
#'     backgroundPosition = "center"
#'   )
#'
#'
#' simplest_df <- data.frame(
#'   a = c(rep("a",9)),
#'   value = c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
#' )
#'
#' # or with a very simple data frame
#' DT::datatable(simplest_df) %>%
#'   formatStyle(
#'     'value',
#'     background = styleColorBar_divergent(simplest_df$value,
#'                                          scales::alpha("forestgreen", 0.4),
#'                                          scales::alpha("gold", 0.4)),
#'     backgroundSize = "100% 90%",
#'     backgroundRepeat = "no-repeat",
#'     backgroundPosition = "center"
#'   )
styleColorBar_divergent <- function(data,
                                    color_pos,
                                    color_neg) {

  max_val <- max(abs(data))
  JS(
    sprintf(
      "isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
      max_val, color_pos, max_val, color_pos, color_neg, color_neg, max_val, max_val))
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
#' a <- 1:9
#' pal <- RColorBrewer::brewer.pal(9,"Set1")
#' map2color(a, pal)
#' plot(a, col = map2color(a, pal), pch = 20, cex = 4)
#'
#' b <- 1:50
#' pal2 <- grDevices::colorRampPalette(
#'   RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50)
#' plot(b, col = map2color(b, pal2), pch = 20, cex = 3)
map2color <- function(x, pal, limits = NULL) {
  if (is.null(limits))
    limits <- range(x)
  pal[findInterval(x, seq(limits[1],
                          limits[2],
                          length.out = length(pal) + 1),
                   all.inside = TRUE)]
}


#' Check colors
#'
#' Check correct specification of colors
#'
#' This is a vectorized version of [grDevices::col2rgb()]
#'
#' @param x A vector of strings specifying colors
#'
#' @return A vector of logical values, one for each specified color - `TRUE` if
#' the color is specified correctly
#' @export
#'
#' @examples
#' # simple case
#' mypal <- c("steelblue", "#FF1100")
#' check_colors(mypal)
#' mypal2 <- rev(
#'   scales::alpha(
#'     colorRampPalette(RColorBrewer::brewer.pal(name = "RdYlBu", 11))(50), 0.4))
#' check_colors(mypal2)
#' # useful with long vectors to check at once if all cols are fine
#' all(check_colors(mypal2))
check_colors <- function(x) {
  sapply(x, function(col) {
    tryCatch(is.matrix(col2rgb(col)),
             error = function(e) FALSE)
  })
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
#' data(res_de_macrophage, package = "GeneTonic")
#' head(res_macrophage_IFNg_vs_naive)
#' res_df <- deseqresult2df(res_macrophage_IFNg_vs_naive)
#' head(res_df)
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

GeneTonic_footer <- fluidRow(
  column(
    width = 1,
    align = "right",
    a(
      href = "https://github.com/federicomarini/GeneTonic",
      target = "_blank",
      img(src = "GeneTonic/GeneTonic.png", height = "50px")
    )
  ),
  column(
    width = 11,
    align = "center",
    "GeneTonic is a project developed by ",
    tags$a(href = "https://federicomarini.github.io", "Federico Marini"),
    " in the Bioinformatics division of the ",
    tags$a(href = "http://www.unimedizin-mainz.de/imbei", "IMBEI"),
    "- Institute for Medical Biostatistics, Epidemiology and Informatics", br(),
    "License: ", tags$a(href = "https://opensource.org/licenses/MIT", "MIT"),
    "- The GeneTonic package is developed and available on ",
    tags$a(href = "https://github.com/federicomarini/GeneTonic", "GitHub")
  )
)

# Shiny resource paths ----------------------------------------------------

.onLoad <- function(libname, pkgname) {
  # Create link to logo
  shiny::addResourcePath("GeneTonic", system.file("www", package = "GeneTonic"))
}


# Some constant values ----------------------------------------------------

.actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC"
.helpbutton_biocstyle <- "color: #0092AC; background-color: #FFFFFF; border-color: #FFFFFF"

gt_downloadButton <- function(outputId, label = "Download", icon = "magic", class = NULL, ...) {
  aTag <- tags$a(id = outputId,
                 class = paste("btn btn-default shiny-download-link", class),
                 href = "",
                 target = "_blank",
                 download = NA,
                 icon(icon),
                 label)
}


.biocgreen <- "#0092AC"
