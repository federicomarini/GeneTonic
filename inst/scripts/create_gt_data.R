library("GeneTonic")
library("macrophage")
data(gse)

# dds object -------------------------------------------------------------------
library("DESeq2")
dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
# no need to save this one, can be readily generated

# annotation object ------------------------------------------------------------
library("org.Hs.eg.db")
anno_df <- data.frame(
  gene_id = rownames(dds_macrophage),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds_macrophage), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds_macrophage)
)

# res object -------------------------------------------------------------------
keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]
library("org.Hs.eg.db")
dds_macrophage <- DESeq(dds_macrophage)
# vst_macrophage <- vst(dds_macrophage)
res_macrophage_IFNg_vs_naive <- results(dds_macrophage,
                                        contrast = c("condition", "IFNg", "naive"),
                                        lfcThreshold = 1, alpha = 0.05)
res_macrophage_IFNg_vs_naive$SYMBOL <- rowData(dds_macrophage)$SYMBOL
library("AnnotationDbi")
# de_symbols_IFNg_vs_naive <- res_macrophage_IFNg_vs_naive[(!(is.na(res_macrophage_IFNg_vs_naive$padj))) & (res_macrophage_IFNg_vs_naive$padj <= 0.05), "SYMBOL"]
de_symbols_IFNg_vs_naive <- deseqresult2df(res_macrophage_IFNg_vs_naive, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds_macrophage)$SYMBOL[rowSums(counts(dds_macrophage)) > 0]
save(res_macrophage_IFNg_vs_naive, file = "data/res_de_macrophage.RData", compress = "xz")

# res_enrich object ------------------------------------------------------------
library("topGO")
topgoDE_macrophage_IFNg_vs_naive <-
  pcaExplorer::topGOtable(de_symbols_IFNg_vs_naive,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "symbol",
                          topTablerows = 500)
write.table(topgoDE_macrophage_IFNg_vs_naive,
            "inst/extdata/topgotable_res_IFNg_vs_naive.txt",
            sep = "\t")
save(topgoDE_macrophage_IFNg_vs_naive, file = "data/res_enrich_macrophage.RData", compress = "xz")
topgoDE_macrophage_IFNg_vs_naive <-
  read.table(system.file("extdata", "topgotable_res_IFNg_vs_naive.txt", package = "GeneTonic"),
             stringsAsFactors = FALSE)


# gostres object ----------------------------------------------------------
library("gprofiler2")
degenes <- deseqresult2df(res_macrophage_IFNg_vs_naive, FDR = 0.01)$SYMBOL
gostres_macrophage <- gost(
  query = degenes, 
  organism = "hsapiens", 
  ordered_query = FALSE, 
  multi_query = FALSE, 
  significant = FALSE, 
  exclude_iea = TRUE, 
  measure_underrepresentation = FALSE, 
  evcodes = TRUE, 
  user_threshold = 0.05, 
  correction_method = "g_SCS", 
  domain_scope = "annotated", 
  # custom_bg =  res_macrophage_IFNg_vs_naive$SYMBOL, 
  numeric_ns = "", 
  sources = "GO:BP", 
  as_short_link = FALSE)
  # as_short_link = TRUE) # alternative, to work on the textual output of the web interface

save(gostres_macrophage, file = "data/gostres_macrophage.RData", compress = "xz")


# fgseaRes object ---------------------------------------------------------
library("dplyr")
library("tibble")
library("fgsea")
res2 <- res_macrophage_IFNg_vs_naive %>%
  as.data.frame() %>% 
  dplyr::select(SYMBOL, stat)
de_ranks <- deframe(res2)
head(de_ranks, 20)
pathways_gmtfile <- gmtPathways("../MSigDBMaker/msigdb_v7.0_files_to_download_locally/msigdb_v7.0_GMTs/c5.bp.v7.0.symbols.gmt")
fgseaRes <- fgsea(pathways = pathways_gmtfile, 
                  stats = de_ranks, 
                  nperm=100000)
fgseaRes <- fgseaRes %>% 
  arrange(desc(NES))
save(fgseaRes, file = "data/fgseaRes.RData", compress = "xz")


# enrichr_output_macrophage -----------------------------------------------
library("enrichR")
dbs <- c("GO_Molecular_Function_2018",
         "GO_Cellular_Component_2018",
         "GO_Biological_Process_2018",
         "KEGG_2019_Human",
         "Reactome_2016",
         "WikiPathways_2019_Human")
degenes <- (deseqresult2df(res_macrophage_IFNg_vs_naive, FDR = 0.01)$SYMBOL)
enrichr_output_macrophage <- enrichr(degenes, dbs)

save(enrichr_output_macrophage, 
     file = "data/enrichr_output_macrophage.RData", 
     compress = "xz")

