
enhance_table <- function(res_enrich,
                          res_de,
                          n_gs = 50,
                          genes_colname = "genes",
                          genesetname_colname = "Term",
                          genesetid_colname = "GO.ID",
                          annotation_obj) {

  # res_enrich has to have a column called containing the genes annotated to the term
  # TODOTODO

  # verify the genesets are sorted in a meaningful way?
  #TODOTODO

  n_gs <- min(n_gs, nrow(res_enrich))

  gs_fulllist <- lapply(seq_len(n_gs), function(go) {
    genes_thisset <- res_enrich[[genes_colname]][go]
    genes_thisset <- unlist(strsplit(genes_thisset,","))

    genesid_thisset <- annotation_obj$gene_id[match(genes_thisset,annotation_obj$gene_name)]

    res_thissubset <- res_de[genesid_thisset,]
    res_thissubset$gene_name <- genes_thisset
    res_thissubset$goterm <- as.factor(res_enrich[[genesetname_colname]][go])
    res_thissubset$gotermshort <- substr(res_enrich[[genesetname_colname]][go],1,50)
    res_thissubset$goid <- res_enrich[[genesetid_colname]][go]
    return(as.data.frame(res_thissubset))
  })
  gs_fulllist <- do.call(rbind, gs_fulllist)

  this_contrast <- (sub(".*p-value: (.*)","\\1",mcols(res_de, use.names=TRUE)["pvalue","description"]))

  # gs_fulllist <- gs_fulllist[nrow(gs_fulllist):1,]
  # gs_fulllist$goterm <- factor(rev(gs_fulllist$goterm))

  p <- ggplot(
    gs_fulllist, aes(
      x = log2FoldChange,
      y = goterm,
      fill = goid,
      text = gene_name
    )) +
    ggtitle(paste0(this_contrast," - TODOTODO")) +

    geom_point(alpha = 0.7, shape = 21, size = 3)  +
    theme_minimal() +
    geom_vline(aes(xintercept = 0), col = "steelblue", alpha = 0.4) +
    theme(legend.position = "none") +
    scale_y_discrete(name = "",
                     labels = substr(as.character(unique(gs_fulllist$goterm)),1,60))

  return(p)
}


# TODOTODO: Z score calculation
