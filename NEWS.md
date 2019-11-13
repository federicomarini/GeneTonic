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
