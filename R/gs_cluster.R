#' Markov Clustering (MCL) for community detection
#'
#' This function implements the Markov Clustering (MCL) algorithm for finding community
#' structure, in an analogous way to other existing algorithms in `igraph`.
#'
#' This implementation has been driven by the nice explanations provided in
#' * https://sites.cs.ucsb.edu/~xyan/classes/CS595D-2009winter/MCL_Presentation2.pdf
#' * https://medium.com/analytics-vidhya/demystifying-markov-clustering-aeb6cdabbfc7
#' * https://github.com/GuyAllard/markov_clustering (python implementation)
#'
#' More info on the MCL: https://micans.org/mcl/index.html, and
#' https://micans.org/mcl/sec_description1.html
#'
#' @references van Dongen, S.M., Graph clustering by flow simulation (2000) PhD thesis,
#' Utrecht University Repository - https://dspace.library.uu.nl/handle/1874/848
#' @references  Enright AJ, van Dongen SM, Ouzounis CA, An efficient algorithm for
#' large-scale detection of protein families (2002) Nucleic Acids Research, Volume 30,
#' Issue 7, 1 April 2002, Pages 1575–1584, https://doi.org/10.1093/nar/30.7.1575
#'
#' @param g The input graph object
#' @param add_self_loops Logical, whether to add self-loops to the matrix by
#' setting the diagonal to `loop_value`
#' @param loop_value Numeric, the value to use for self-loops
#' @param mcl_expansion Numeric, cluster expansion factor for the Markov clustering
#' iteration - defaults to 2
#' @param mcl_inflation Numeric, cluster inflation factor for the Markov clustering
#' iteration - defaults to 2
#' @param allow_singletons Logical; if `TRUE`, single isolated vertices are allowed
#' to form their own cluster. If set to `FALSE`, all clusters of size = 1 are
#' grouped in one cluster (to be interpreted as background noise).
#' @param max_iter Numeric value for the maximum number of iterations for the
#' Markov clustering
#' @param return_node_names Logical, if the graph is named and set to `TRUE`, returns
#' the node names.
#' @param return_esm Logical, controlling whether the equilibrium state matrix should be returned
#'
#' @return This function returns a `communities` object, containing the numbers of
#' the assigned membership (in the slot `membership`). Please see the
#' [igraph::communities()] manual page for additional details
#'
#' @export
#'
#' @examples
#' library("igraph")
#' g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
#' g <- add_edges(g, c(1, 6, 1, 11, 6, 11))
#' cluster_markov(g)
#' V(g)$color <- cluster_markov(g)$membership
#' plot(g)
cluster_markov <- function(g,
                           add_self_loops = TRUE,
                           loop_value = 1,
                           mcl_expansion = 2,
                           mcl_inflation = 2,
                           allow_singletons = TRUE,
                           max_iter = 100,
                           return_node_names = TRUE,
                           return_esm = FALSE) {
  # g must be a graph
  if (!is(g, "igraph")) {
    stop("You need to provide an igraph object as input")
  }

  stopifnot(is.logical(add_self_loops))
  stopifnot(loop_value >= 0)
  stopifnot(mcl_expansion > 1)
  stopifnot(mcl_inflation > 1)
  stopifnot(loop_value >= 0)
  stopifnot(max_iter > 0)

  # graph to adjacency matrix
  adj_mat <- igraph::as_adjacency_matrix(g)
  adj_mat <- as.matrix(adj_mat) # to enforce non-sparse matrix

  converged <- FALSE

  # Initialize self-loops
  if (add_self_loops) {
    diag(adj_mat) <- loop_value
  }

  # Normalize (sum by col must be 1)
  adj_mat_norm <- t(adj_mat / colSums(adj_mat))

  cur_mat_norm <- adj_mat_norm
  for (i_mcl in seq_len(max_iter)) {
    # message(i_mcl)
    last_mat <- cur_mat_norm

    exp_mat <- cur_mat_norm %^% mcl_expansion
    inf_mat <- exp_mat^mcl_inflation
    # inf_mat_norm <- t(inf_mat/colSums(inf_mat))
    inf_mat_norm <- apply(inf_mat, MARGIN = 2, FUN = function(matcol) {
      matcol / sum(matcol)
    })

    ## TODO: optional pruning?
    # idea: inspect matrix and set small values directly to zero (assume they would have reached there eventually anyways).
    if (identical(inf_mat_norm, last_mat)) {
      converged <- TRUE
      break
    }

    cur_mat_norm <- inf_mat_norm
  }

  if (converged & is.na(cur_mat_norm[1, 1])) {
    stop("An error occurred after convergence - maybe you set `add_self_loops` to FALSE?")
  }

  # getting the attractors - non-zero elements of the matrix diagonal
  clu_attractors <- which(diag(cur_mat_norm) > 0)
  # store the clusters
  clusters <- vector(mode = "list", length(clu_attractors))

  # the nodes in the same row as each attractor form a cluster
  for (i in seq_along(clu_attractors)) {
    cur_att <- clu_attractors[i]
    cur_members <- which(cur_mat_norm[cur_att, ] > 0)
    clusters[[i]] <- cur_members
  }

  # chop off the identical ones
  clusters <- unique(clusters)

  # from clusters sets to membership as label
  clu_memberships <- rep(NA, nrow(adj_mat))
  for (i in seq_along(clusters)) {
    this_cluster <- clusters[[i]]
    clu_memberships[this_cluster] <- i
  }

  # handling the singletons
  if (!allow_singletons) {
    dub <- duplicated(clu_memberships) + duplicated(clu_memberships, fromLast = TRUE)
    n_singletons <- sum(dub == 0)
    clu_memberships[dub == 0] <- 0
    # reshifting to not having gaps in the cluster numbers
    clu_memberships[dub != 0] <- clu_memberships[dub != 0] - n_singletons
  }

  res <- list()
  if (return_node_names & igraph::is_named(g)) {
    res$names <- V(g)$name
  }
  res$vcount <- vcount(g)
  res$algorithm <- "mcl"
  res$iterations <- i_mcl - 1
  res$membership <- clu_memberships
  if (return_esm) {
    res$esm <- cur_mat_norm
  }

  class(res) <- "communities" # to take advantage of the goodies for printing and co

  return(res)
}


# nice to test on the set of igraph examples - say, louvain
