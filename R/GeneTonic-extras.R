#' Title
#'
#' TODO
#'
#' @param x TODO
#' @param y TODO
#'
#' @return TODO
#' @export
#'
#' @seealso https://en.wikipedia.org/wiki/Overlap_coefficient
#'
#' @examples
#' # TODO
overlap_coefficient <- function (x, y)
{
  length(intersect(x, y))/min(length(x), length(y))
}

#' Title
#'
#' TODO
#'
#' @param x TODO
#' @param y TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
overlap_jaccard_index <- function (x, y)
{
  length(intersect(x, y))/length(unique(c(x, y)))
  # about 2x faster than using union()
}



#' Title
#'
#' TODO
#'
#' @param x TODO
#' @param pal TODO
#' @param limits TODO
#'
#' @return TODO
#' @export
#'
#' @examples
#' # TODO
map2color <- function(x, pal, limits = NULL) {
  if(is.null(limits))
    limits=range(x)
  pal[findInterval(x, seq(limits[1],
                          limits[2],
                          length.out=length(pal)+1),
                   all.inside=TRUE)]
}
