`GeneTonic` requires you to provide as input a `GeneTonicList`, i.e. a list-like object, storing four essential components, based on the widely used `DESeq2` framework - if you are using alternative workflows, you might need to convert them beforehand.  
A `GeneTonicList` contains these named slots:

1. `dds`, a `DESeqDataSet` object stoging all the information related to the expression matrix.
2. `res_de`, the `DESeqResults` object with the results of Differential Expression analysis.
3. `res_enrich`, a simple data frame with the results from functional enrichment analysis tools (see the documentation to find out which tools are supported via the `shaker_` functions). Column names must follow some specific requirements to guarantee interoperability within `GeneTonic`.
4. `annotation_obj`, the annotation data frame, composed at least of two columns, `gene_id`, with a set of unambiguous identifiers (e.g. ENSEMBL ids) corresponding to the row names of the `dds` object, and `gene_name`, containing e.g. HGNC-based gene symbols. 

You can create a `GeneTonicList` by calling

```
mygtlobject <- GeneTonicList(
  dds = dds,
  res_de = res_de,
  res_enrich = res_enrich,
  annotation_obj = anno_df
)
```

If desired, you can export the `GeneTonicList` to a serialized object, and then upload this in the current panel of the `GeneTonic` app:

```
saveRDS(mygtlobject, file = "mygtlobject.RDS")
```

If running `GeneTonic` from the command line, you can simply execute these lines to run a full example on the `macrophage` dataset:

```
library("GeneTonic")
example("GeneTonic", ask = FALSE)
```


<hr>

For additional information about the format, you can refer to

* the `GeneTonic` package vignette (to be found at http://bioconductor.org/packages/GeneTonic/)
* the original `GeneTonic` publication (https://doi.org/10.1186/s12859-021-04461-5, on BMC Bioinformatics), in particular the workflow depicted in Figure 1
* a series of protocols to appear in *Current Protocols in Bioinformatics*, where `GeneTonic` is featured, together with the `pcaExplorer` and `ideal` packages

<div align="center">
<img src="GeneTonic/GeneTonic.png" alt="" width="150" />
</div>
