#' Title TODO
#'
#' TODO
#'
#' @param res_enrich TODO
#' @param res_de TODO
#' @param annotation_obj TODO
#' @param n_gs TODO
#' @param genes_colname TODO
#' @param genesetname_colname TODO
#' @param genesetid_colname TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
gs_alluvial <- function(res_enrich,
                        res_de,
                        annotation_obj,
                        n_gs = 5,
                        genes_colname = "genes",
                        genesetname_colname = "Term",
                        genesetid_colname = "GO.ID") {

  # res_enhanced <- get_aggrscores(res_enrich, res_de, annotation_obj = annotation_obj)

  enriched_gsids <- res_enrich[[genesetid_colname]]
  enriched_gsnames <- res_enrich[[genesetname_colname]]
  enriched_gsdescs <- vapply(enriched_gsids,
                             function(arg) Definition(GOTERM[[arg]]),
                             character(1))

  enrich2list <- lapply(seq_len(n_gs), function(gs) {
    # goterm <- res_enrich$Term[gs]
    go_genes <- res_enrich$genes[gs]
    go_genes <- strsplit(go_genes, ",") %>% unlist
    return(go_genes)
  })
  names(enrich2list) <- enriched_gsnames[seq_len(n_gs)]

  list2df <- lapply(seq_len(length(enrich2list)), function(gs) {
    data.frame(
      gsid = rep(names(enrich2list[gs]), length(enrich2list[[gs]])),
      gene = enrich2list[[gs]])
  })
  list2df <- do.call("rbind", list2df)

  list2df$value <- 1
  colnames(list2df) <- c("source", "target", "value")
  list2df$source <- as.character(list2df$source)
  list2df$target <- as.character(list2df$target)
  list2df$target <- paste(list2df$target, " ", sep = "") # genes might need a space for better rendering...

  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
    name = c(as.character(list2df$source), as.character(list2df$target)) %>% unique()
  )

  # connections must be provided using id, not using real name like in the links data.frame...
  list2df$IDsource <- match(list2df$source, nodes$name) - 1
  list2df$IDtarget <- match(list2df$target, nodes$name) - 1

  # list2df %>% head

  allnodes <- c(unique(list2df$source), unique(list2df$target))
  allcols <- c(
    viridis(length(unique(list2df$source))),           # for genesets
    rep("steelblue", length(unique(list2df$target)))    # for genes
  )

  p <- plot_ly(
    type = "sankey",
    domain = list(
      x =  c(0,1),
      y =  c(0,1)
    ),
    orientation = "h",
    valueformat = ".0f",
    # valuesuffix = " genes in the set",  # TODO: should be done conditional on what the node is

    node = list(
      label = allnodes,
      color = allcols,
      pad = 15,
      thickness = 20,
      line = list(color = "black", width = 0.5)
    ),

    link = list(
      source = list2df$IDsource,
      target = list2df$IDtarget,
      value =  list2df$value
    )
  ) %>%
    plotly::layout(
      title = "Geneset-Gene Sankey Diagram", font = list(size = 10)
    )
  return(p)
}

