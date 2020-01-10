# GeneTonic 0.9.0

## New features

* `GeneTonic` sports a blazing new hex sticker - say bye to the original draft!
* The overview DT `datatable`s has some styling with color bars - e.g. for DE results - to enhance the visual perception of numeric values (e.g. log2FoldChange)
* `gs_heatmap` can now take a custom list of gene identifiers (when no geneset is passed)
* The color palettes in enrichment maps now respect the values and the range specified of the numeric values to be used for mapping to colors
* `gs_mds` is now optionally returning a data.frame, to be further used for custom plotting or downstream processing
* `gs_summary_overview` now has coloring enabled by the variable of choice

## Other notes

* The UI has received some restyling (e.g. in the choice of the icons for the dropdown menus, or the name of some buttons)
* Added tour contents for most of the functionality
* Added link to the demo instance
* Added examples for overlap functions, gene info buttons, map2color, and deseqresult2df
* Extended documentation of some parameters
* Some functions have gained an alias for calling them: `gs_spider` is equivalent to `gs_radar`, and `gs_sankey` is equivalent to `gs_alliuvial`

# GeneTonic 0.8.0

## New features

* `GeneTonic` now delivers bundled example objects to make examples and tests slim
* `gs_volcano` can now plot points by different colors according to the columns of interest
* `GeneTonic` has a fully fledged manual describing its functionality and user interface

## Other notes

* Now using ids for genes and genesets for exchanging information in the app
* Added examples for all functions
* Most tabs have working tours - anchor and text elements

# GeneTonic 0.7.0

## New features

* Introduced a uniform interface for calculating different similarity/distance matrices. This enables the usage in the different functions that might need such matrices for further downstream processing (e.g. `enrichment_map()`, `gs_mds()`)
* First appearance of `gs_dendro()` to display distance matrices with some visualization sugar, as an alternative to other methods
* The `n_gs` and `gs_ids` are exposed to more functions to enable custom subsets of the enrichment results to be inspected

# GeneTonic 0.6.0

## New features

* `gs_heatmap` now relies on `ComplexHeatmap`, to avoid the issues with Shiny of not displaying the outputs in the app, and enabling a comfortable heatmap annotation
* Many functions gain the possibility to pass a set of custom geneset identifiers to be added to the top N sets (default): among these, `gs_mds`, `gs_volcano` (parameter: `gs_labels`), `gs_alluvial`, `ggs_network`, `enrichment_map`, and `enhance_table` (using `gs_ids`)

## Other notes

* `gs_ggheatmap` got renamed to `gs_scoresheat`
* The report generated from the bookmarked content is expanded in its default content

# GeneTonic 0.5.0

## New features

* `GeneTonic` now enforces a format for `res_enrich`, and provides some conversion functions, `shake_*()`. Requirements are specified in the documentation, if an appropriate converter does not (yet) exist.
* The reporting feature is active to some extent on the bookmarked content.

# GeneTonic 0.4.0

## New features

* Added functionality for bookmarking
* Bookmarking can work (PoP) by pressing a key (left control)!
* `gene_plot` can enforce a plot type overriding the default based on the number of samples per condition
* `GeneTonic` uses now `bs4Dash` and many of its nice features, replacing the previous implementation based on `shinydashboard`

# GeneTonic 0.3.0

## Other notes

* Rearranging the order of parameters to harmonize it across functions, and uniforming similarly called parameters

# GeneTonic 0.2.0

## New features

* added geneset scoring calculation and corresponding heatmap

# GeneTonic 0.1.0

## New features

* much of the functionality available, in a proof of concept format

# GeneTonic 0.0.1

## New features

* backbone of the project started!
