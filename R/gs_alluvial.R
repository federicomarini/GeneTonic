#' Title
#'
#' @param res_enrich
#' @param res_de
#' @param annotation_obj
#' @param n_gs
#' @param genes_colname
#' @param genesetname_colname
#' @param genesetid_colname
#' @param ggstatic
#'
#' @return
#' @export
#'
#' @examples
gs_alluvial <- function(res_enrich,
                        res_de,
                        annotation_obj,
                        n_gs = 5,
                        genes_colname = "genes",
                        genesetname_colname = "Term",
                        genesetid_colname = "GO.ID",
                        ggstatic = FALSE) {

  res_enhanced <- get_aggrscores(res_enrich, res_de, annotation_obj = annotation_obj)

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

  mygsdf <- list2df

  if (ggstatic) {
    mygsdf$value <- 1

    ggplot(as.data.frame(mygsdf),
           aes(y = value, axis1 = gsid, axis2 = gene)) +
      geom_alluvium(aes(fill = gsid), width = 1/12) +
      geom_stratum(width = 1/12, fill = "grey", color = "black") +
      geom_label(stat = "stratum", label.strata = TRUE) +
      # geom_label_repel(stat = "stratum", label.strata = TRUE) +
      scale_x_discrete(limits = c("Gene set", "Gene"), expand = c(0.5, 0.2)) +
      theme_void() -> p
    return(p)
  } else {

    mygsdf2 <- list2df
    mygsdf2$value <- 1
    colnames(mygsdf2) <- c("source", "target", "value")
    mygsdf2$source <- as.character(mygsdf2$source)
    mygsdf2$target <- as.character(mygsdf2$target)
    mygsdf2$target <- paste(mygsdf2$target, " ", sep="")

    # From these flows we need to create a node data frame: it lists every entities involved in the flow
    nodes <- data.frame(name=c(as.character(mygsdf2$source), as.character(mygsdf2$target)) %>% unique())

    # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
    mygsdf2$IDsource=match(mygsdf2$source, nodes$name)-1
    mygsdf2$IDtarget=match(mygsdf2$target, nodes$name)-1

    mygsdf2 %>% head

    allnodes <- c(unique(mygsdf2$source), unique(mygsdf2$target))

    library(plotly)

    p <- plot_ly(
      type = "sankey",
      domain = list(
        x =  c(0,1),
        y =  c(0,1)
      ),
      orientation = "h",
      valueformat = ".0f",
      valuesuffix = " genes in the set",

      node = list(
        label = allnodes,
        color = c(viridis(length(unique(mygsdf2$source))), rep("steelblue",length(unique(mygsdf2$target)))),
        pad = 15,
        thickness = 20,
        line = list(
          color = "black",
          width = 0.5
        )
      ),

      link = list(
        source = mygsdf2$IDsource,
        target = mygsdf2$IDtarget,
        value =  mygsdf2$value
      )
    ) %>%
      plotly::layout(
        title = "Geneset-Gene Sankey Diagram",
        font = list(
          size = 10
        )
      )
    return(p)
  }

}

