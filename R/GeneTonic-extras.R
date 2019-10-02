#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
overlap_ratio <- function (x, y)
{
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x, y)))
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
