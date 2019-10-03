#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @seealso https://en.wikipedia.org/wiki/Overlap_coefficient
#'
#' @examples
overlap_coefficient <- function (x, y)
{
  length(intersect(x, y))/min(length(x), length(y))
}

overlap_jaccard_index <- function (x, y)
{
  length(intersect(x, y))/length(unique(c(x, y)))
  # about 2x faster than using union()
}



#' Title
#'
#' @param x
#' @param pal
#' @param limits
#'
#' @return
#' @export
#'
#' @examples
map2color <- function(x, pal, limits = NULL) {
  if(is.null(limits))
    limits=range(x)
  pal[findInterval(x, seq(limits[1],
                          limits[2],
                          length.out=length(pal)+1),
                   all.inside=TRUE)]
}
