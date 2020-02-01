## What is `GeneTonic`?

`GeneTonic` is a Bioconductor package whose aim is to analyze and integrate the results from Differential Expression analysis and functional enrichment analysis, together with the original expression matrix.

As a distinctive trait, this app combines the ease and beauty of interactive analysis with the technical robustness and the practicality of generating a reproducible report.

## What do I need to use `GeneTonic`?

Four main ingredients are required:

1. `dds`: a `DESeqDataSet` object, the main component in the `DESeq2` framework, which extends the widely adopted `SummarizedExperiment` class. This object will store the information related to the expression matrix.

2. `res_de`: a `DESeqResults` object, the results of Differential Expression analysis, computed on the `dds` object above. This extends the `S4Vectors::DataFrame` class.

3. `res_enrich`: a `data.frame`, containing the results of the enrichment analysis, generally computed on the basis of the DE results and the expression matrix. The `GeneTonic` main functions require some specific columns to be present, specified in the package documentation. Common formats (e.g. from `pcaExplorer::topGOtable` or from `ClusterProfiler`) are supported with conversion functions.

4. `annotation_obj`: the annotation `data.frame`, composed at least of two columns, `gene_id`, with a set of unambiguous identifiers (e.g. ENSEMBL ids) corresponding to the row names of the `dds` object, and `gene_name`, containing e.g. HGNC-based gene symbols. 

For more detail on these and more aspects (getting to know the user interface, using the functions in a normal R session), please consult the vignette, and try out the tours for each tab in a running instance of the app (e.g. the demo available at http://shiny.imbei.uni-mainz.de:3838/GeneTonic)
