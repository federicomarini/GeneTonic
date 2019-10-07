library(GeneTonic)

library(macrophage)
library(magrittr)
dir <- system.file("extdata", package = "macrophage")
coldata_macrophage <- read.csv(file.path(system.file("extdata", package = "macrophage"), "coldata.csv"))
coldata_macrophage$files <- file.path(system.file("extdata", package = "macrophage"), "quants", coldata_macrophage$names, "quant.sf.gz")

coldata_macrophage$condition <- coldata_macrophage$condition_name
coldata_macrophage$line <- coldata_macrophage$line_id
coldata_macrophage$condition <- relevel(coldata_macrophage$condition, "naive")

head(coldata_macrophage)

library(SummarizedExperiment)
library(tximeta)
se_macrophage <- tximeta(coldata_macrophage, dropInfReps = TRUE)
se_macrophage

# saveRDS(se_macrophage,file="WIP/se_macrophage.rds")
# dir.create("WIP")
# se_macrophage <- readRDS(file = "WIP/se_macrophage.rds")
gse_macrophage <- summarizeToGene(se_macrophage)

assayNames(se_macrophage)
gse_macrophage <- summarizeToGene(se_macrophage)
gse_macrophage
# adding gene names to facilitate readout later
library(org.Hs.eg.db)
gse_macrophage <- addIds(gse_macrophage, "SYMBOL")

library(DESeq2)
dds_macrophage <- DESeqDataSet(gse_macrophage, design = ~line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)

anno_df <- pcaExplorer::get_annotation_orgdb(dds_macrophage, "org.Hs.eg.db", "ENSEMBL")
## using counts and average transcript lengths from tximeta
keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
dds_macrophage <- dds_macrophage[keep, ]

library("org.Hs.eg.db")
dds_macrophage <- addIds(dds_macrophage, "SYMBOL")
dds_macrophage <- DESeq(dds_macrophage)
vst_macrophage <- vst(dds_macrophage)
res_macrophage_IFNg_vs_naive <- results(dds_macrophage,
                                        contrast = c("condition", "IFNg", "naive"),
                                        lfcThreshold = 1, alpha = 0.01)
summary(res_macrophage_IFNg_vs_naive)
res_macrophage_IFNg_vs_naive$SYMBOL <- rowData(dds_macrophage)$SYMBOL

library("AnnotationDbi")

de_symbols_IFNg_vs_naive <- res_macrophage_IFNg_vs_naive[ (!(is.na(res_macrophage_IFNg_vs_naive$padj))) & (res_macrophage_IFNg_vs_naive$padj <= 0.05), "SYMBOL"]
bg_ids <- rowData(dds_macrophage)$SYMBOL[rowSums(counts(dds_macrophage)) > 0]

library("topGO")
topgoDE_macrophage_IFNg_vs_naive <- pcaExplorer::topGOtable(de_symbols_IFNg_vs_naive,
                                                            bg_ids,
                                                            ontology = "BP",
                                                            mapping = "org.Hs.eg.db",
                                                            geneID = "symbol")
