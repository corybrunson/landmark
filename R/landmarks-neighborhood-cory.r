#' @title Neighborhood-based Landmark Sets
#' @author Jason Cory Brunson
#' @description Compute landmark sets based on nearest neighborhoods.
#' @details These functions adapt the maxmin procedure to produce landmark
#'   points dispersed according to the orders in which they are reached from
#'   each other, rather than to their distances from each other. (Say more.)
#'
#'   One, both, or neither of `number` and `cardinality` may be passed values.
#'   If neither is specified, then `cardinality` is internally set to `1L` so
#'   that a complete landmark set is generated. If the values yield
#'   neighborhoods that do not cover `x`, then, effectively, `number` is
#'   increased until the cardinality necessary to cover `x` is at most
#'   `cardinality`.
#' @name landmarks_neighborhood
#' @param x a data matrix.
#' @param y a data matrix, or `NULL` to recycle `x`.
#' @param number the number of landmarks requested.
#' @param cardinality the desired cardinality of each cover set.
#' @param dist_method the distance metric to use. Any distance measure in the
#'   \code{proxy} package is supported.
#' @param seed_index the first landmark to seed the algorithm.
#' @param shuffle_data whether to first randomly shuffle the data.
NULL

#' @rdname landmarks_neighborhood
#' @export
chebycenter_firstlast_cory <- function(
  x, y = NULL,
  dist_method = "euclidean"
) {
  # validate inputs
  stopifnot(is.matrix(x))
  if (is.null(y)) y <- x
  stopifnot(is.matrix(y))

  # obtain the firstlast subset
  # points in `x` that maximize the maximum-rank distance points in `y`
  fl_idx <- seq(nrow(x))
  for (i in rev(seq(nrow(y)))) {
    rk_i <- apply(x[fl_idx, , drop = FALSE], 1, function(r) {
      length(which(rank(proxy::dist(t(r), y, method = dist_method),
                        ties.method = "max") == i))
    })
    fl_idx <- fl_idx[which(rk_i == min(rk_i))]
    if (length(fl_idx) == 1L) break
  }

  fl_idx[[1L]]
}

#' @rdname landmarks_neighborhood
#' @export
landmarks_lastfirst_cory <- function(
  x, number = NULL, cardinality = NULL,
  dist_method = "euclidean", seed_index = 1L, shuffle_data = FALSE
) {
  # validate inputs
  stopifnot(is.matrix(x))
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
  # strike seed and (other) duplicates index from free indices
  perm_idx <- c(lf_idx, free_idx[-lf_idx])
  free_idx[perm_idx[duplicated(x[perm_idx])]] <- 0L
  landmark_rank <- matrix(NA, nrow = nrow(x), ncol = 0)

  # require a number of neighborhoods or a neighborhood cardinality (or both)
  #if (is.null(number) && is.null(cardinality)) number <- length(free_idx)
  if (is.null(number) && is.null(cardinality)) cardinality <- 1L

  # recursively construct landmark set
  for (i in seq_along(free_idx)) {

    # update vector of landmark points
    landmark_idx[[i]] <- lf_idx[[1L]]

    # update vector of free points
    if (free_idx[[landmark_idx[[i]]]] == 0L)
      stop("Landmark choice is a duplicate point.")
    free_idx[[landmark_idx[[i]]]] <- 0L

    # augment ranks from new landmark point
    landmark_rank <- cbind(
      landmark_rank,
      rank(proxy::dist(x,
                       x[landmark_idx[[i]], , drop = FALSE],
                       method = dist_method),
           ties.method = "max")
    )
    # sort the points' rankings
    landmark_rank[] <- t(apply(landmark_rank, 1L, sort))

    # refresh the minimum cardinality
    min_cardinality <- max(landmark_rank[c(free_idx, landmark_idx), 1L])
    if (min_cardinality < ncol(landmark_rank))
      landmark_rank <- landmark_rank[, seq(min_cardinality), drop = FALSE]

    # exhaustion breaks
    if (all(free_idx == 0L)) break
    # parameter breaks
    if (! is.null(number)) {
      # discontinue if the desired number of sets has been reached
      if (i >= number) {
        # discontinue if there is no cardinality requirement
        if (is.null(cardinality)) {
          break
        } else {
          # continue if desired cardinality requires a greater number of sets
          if (min_cardinality <= cardinality) break
        }
      }
    } else {
      # continue if desired cardinality requires more sets
      if (min_cardinality <= cardinality) break
    }

    # obtain the lastfirst subset
    lf_idx <- free_idx[free_idx != 0L]
    for (j in seq(ncol(landmark_rank))) {
      lf_idx <- lf_idx[landmark_rank[lf_idx, j] ==
                         max(landmark_rank[lf_idx, j])]
      if (length(lf_idx) == 1L) break
    }

  }

  landmark_idx[seq(i)]
}
