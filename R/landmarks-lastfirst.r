#' @title Neighborhood-based Landmark Sets
#' @author Jason Cory Brunson
#' @author Yara Skaf
#' @description Compute landmark sets based on nearest neighborhoods.
#' @details These functions adapt the maxmin procedure to produce landmark
#'   points dispersed according to the orders in which they are reached from
#'   each other, rather than to their distances from each other. (Say more.)
#'
#'   One, both, or neither of `num_sets` and `cardinality` may be passed values.
#'   If neither is specified, then `cardinality` is internally set to `1L` so
#'   that a complete landmark set is generated. If the values yield
#'   neighborhoods that do not cover `x`, then, effectively, `num_sets` is
#'   increased until the cardinality necessary to cover `x` is at most
#'   `cardinality`.
#' @name landmarks_lastfirst
#' @param x a data matrix.
#' @param dist_method a character string specifying the distance metric to use;
#'   passed to `proxy::dist(method)`. Any distance measure in the \code{proxy}
#'   package is supported.
#' @param ties_method a character string specifying the method for handling
#'   ties; passed to `rank(ties.method)`. Only `"min"` has been tested and is
#'   recommended.
#' @param pick_method a character string specifying the method for selecting
#'   from indistinguishable points, either `"first"` (the default) or
#'   `"random"`.
#' @param num_sets a positive integer; the desired number of landmark points, or
#'   of sets in a neighborhood cover.
#' @param cardinality a positive integer; the desired cardinality of each
#'   landmark neighborhood, or of each set in a landmark cover.
#' @param seed_index an integer (the first landmark to seed the algorithm) or
#'   one of the character strings `"random"` (to select a seed uniformly at
#'   random) and `"firstlast"` (to select a seed from the firstlast set).
#' @param shuffle_data logical; whether to randomly shuffle the data before
#'   running a procedure.
NULL

#' @rdname landmarks_lastfirst
#' @export
firstlast_R <- function(
  x,
  dist_method = "euclidean", ties_method = "min"
) {

  # initialize colex-minimum out-rank-distance sequence
  seq_min <- rep(nrow(x), nrow(x))
  # initialize firstlast set
  fl_idx <- integer(0)

  # across all points
  for (idx in seq(nrow(x))) {

    # out-rank-distance sequence
    seq_idx <- sort(rank(proxy::dist(x[idx, , drop = FALSE], x,
                                     method = dist_method),
                         ties.method = ties_method))
    # latest rank at which it disagrees with the reigning minimum sequence
    diff_last <- suppressWarnings(max(which(seq_idx != seq_min)))

    if (diff_last == -Inf) {
      # if equal to reigning minimum sequence, append to firstlast set
      fl_idx <- c(fl_idx, idx)
    } else if (seq_idx[[diff_last]] < seq_min[[diff_last]]) {
      # if less than reigning minimum sequence, replace and reinitialize
      seq_min <- seq_idx
      fl_idx <- c(idx)
    }

  }

  # return firstlast subset
  fl_idx
}

#' @rdname landmarks_lastfirst
#' @export
landmarks_lastfirst_R <- function(
  x,
  dist_method = "euclidean", ties_method = "min", pick_method = "first",
  num_sets = NULL, cardinality = NULL,
  seed_index = 1L, shuffle_data = FALSE
) {
  # validate inputs
  stopifnot(is.matrix(x))
  stopifnot(toupper(dist_method) %in% toupper(proxy::pr_DB$get_entry_names()))

  # handle seed selection
  if (is.character(seed_index)) {
    seed_index <- switch (
      match.arg(seed_index, c("random", "firstlast")),
      random = sample(nrow(x), size = 1L),
      firstlast = firstlast_R(x,
                              dist_method = dist_method,
                              ties_method = ties_method)
    )
  }
  stopifnot(seed_index >= 1L && seed_index <= nrow(x))

  # initialize lastfirst, free, and landmark index sets
  lf_idx <- seed_index
  lmk_idx <- vector(mode = "integer", nrow(x))
  free_idx <- seq(nrow(x))
  # strike seed and (other) duplicates index from free indices
  perm_idx <- c(lf_idx, free_idx[-lf_idx])
  free_idx[perm_idx[duplicated(x[perm_idx])]] <- 0L
  lmk_rank <- matrix(NA, nrow = nrow(x), ncol = 0)

  # require a number of neighborhoods or a neighborhood cardinality (or both)
  #if (is.null(num_sets) && is.null(cardinality)) num_sets <- length(free_idx)
  if (is.null(num_sets) && is.null(cardinality)) cardinality <- 1L

  # recursively construct landmark set
  for (i in seq_along(free_idx)) {

    # update vector of landmark points
    lmk_idx[[i]] <- lf_idx[[switch(
      match.arg(pick_method, c("first", "random")),
      first = 1L,
      random = sample(length(lf_idx), 1L)
    )]]

    # update vector of free points
    if (free_idx[[lmk_idx[[i]]]] == 0L)
      stop("A duplicate landmark point was selected, in error.")
    free_idx[[lmk_idx[[i]]]] <- 0L

    # augment in-ranks from new landmark point
    lmk_rank <- cbind(
      # each row contains the in-ranks from previous landmark points
      lmk_rank,
      # in-ranks of all points from newest landmark point
      rank(proxy::dist(x[lmk_idx[[i]], , drop = FALSE], x,
                       method = dist_method),
           ties.method = ties_method)
    )
    # sort each available point's in-ranks to the landmark points
    lmk_rank[] <- t(apply(lmk_rank, 1L, sort))

    # refresh the minimum cardinality
    min_card <- max(lmk_rank[c(free_idx, lmk_idx), 1L])
    # drop ranks past minimum cardinality
    if (min_card < ncol(lmk_rank))
      lmk_rank <- lmk_rank[, seq(min_card), drop = FALSE]

    # exhaustion breaks
    if (all(free_idx == 0L)) break
    # parameter breaks
    if (! is.null(num_sets)) {
      # discontinue if the desired number of sets has been reached
      if (i >= num_sets) {
        # discontinue if there is no cardinality requirement
        if (is.null(cardinality)) {
          break
        } else {
          # continue if desired cardinality requires a greater number of sets
          if (min_card <= cardinality) break
        }
      }
    } else {
      # continue if desired cardinality requires more sets
      if (min_card <= cardinality) break
    }

    # obtain the lastfirst subset
    lf_idx <- free_idx[free_idx != 0L]
    for (j in seq(ncol(lmk_rank))) {
      # points with maximum revlex rank-in-distance counts
      # = points with minimum lex rank-in-distance counts
      # = points with maximum lex rank-in-distance sequence
      lf_idx <- lf_idx[lmk_rank[lf_idx, j] == max(lmk_rank[lf_idx, j])]
      if (length(lf_idx) == 1L) break
    }

  }

  # print warnings if a parameter was adjusted
  if (! is.null(num_sets)) {
    if (i > num_sets) {
      warning("Cover required ", i, " (> num_sets = ", num_sets, ") ",
              "sets of cardinality ", cardinality, ".")
    } else if (i < num_sets) {
      warning("Only ", i, " (< num_sets = ", num_sets, ") ",
              "distinct landmark points were found.")
    }
  }

  lmk_idx[seq(i)]
}
