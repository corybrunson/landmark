#' @title Firstlast-generated Landmarks
#' @author Jason Cory Brunson
#' @description Computes a landmark set based on nearest neighborhoods.
#' @details This function adapts the maxmin procedure to produce landmark points
#'   dispersed according to the orders in which they are reached from each
#'   other, rather than to their distances from each other. (Say more.)
#' @name landmarks_lastfirst
#' @param x a data matrix.
#' @param number the number of landmarks requested.
#' @param cardinality the desired cardinality of each cover set.
#' @param dist_method the distance metric to use. Any distance measure in the
#'   \code{proxy} package is supported.
#' @param seed_index the first landmark to seed the algorithm.
#' @param shuffle_data whether to first randomly shuffle the data.
#' @export
landmarks_lastfirst_cory <- function(
  x, number = NULL, cardinality = NULL,
  dist_method = "euclidean", seed_index = 1L, shuffle_data = FALSE
) {
  # validate inputs
  stopifnot(is.matrix(x))
  # require a number of neighborhoods or a neighborhood cardinality (or both)
  stopifnot(! is.null(number) || ! is.null(cardinality))
  stopifnot(toupper(dist_method) %in% toupper(proxy::pr_DB$get_entry_names()))
  stopifnot(seed_index >= 1L && seed_index <= nrow(x))

  # shuffle data if desired
  if (shuffle_data){
    x <- x[sample(nrow(x)), , drop = FALSE]
  }

  # initialize lastfirst, free, and landmark index sets
  lf_idx <- seed_index
  landmark_idx <- vector(mode = "integer", nrow(x))
  free_idx <- seq(nrow(x))
  landmark_rank <- matrix(NA, nrow = nrow(x), ncol = 0)

  # recursively construct landmark set
  for (i in seq(nrow(x))) {

  }

  landmark_idx[seq(i)]
}
