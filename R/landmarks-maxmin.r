#' @title Ball-based Landmark Sets
#' @author Matt Piekenbrock
#' @author Jason Cory Brunson
#' @author Yara Skaf
#' @description Compute landmark sets based on fixed-radius balls.
#' @details This function uses the maxmin procedure to produce a set of evenly
#'   spaced landmark points from a data set. Maxmin is a simple greedy algorithm
#'   that is relatively efficient, but it has a tendency to pick out extremal
#'   points.
#'
#'   One, both, or neither of `num_sets` and `radius` may be passed values. If
#'   neither is specified, then `num_sets` is defaulted to `24L`. If the values
#'   yield balls that do not cover `x`, then their number is increased until the
#'   radius necessary to cover `x` is at most `radius`. To generte a complete
#'   landmark set, use `radius = 0L`.
#' @name landmarks_maxmin
#' @param x a data matrix.
#' @param dist_method a character string specifying the distance metric to use;
#'   passed to `proxy::dist(method)`. Any distance measure in the \code{proxy}
#'   package is supported.
#' @param pick_method a character string specifying the method for selecting one
#'   among indistinguishable points, either `"first"` (the default), `"last"`,
#'   or `"random"`.
#' @param num_sets a positive integer; the desired number of landmark points, or
#'   of sets in a ball cover.
#' @param radius a positive number; the desired radius of each
#'   landmark ball, or of each set in a ball cover.
#' @param frac logical; whether to treat `radius` as a fraction of the diameter
#'   of `x`.
#' @param seed_index an integer (the first landmark to seed the algorithm) or
#'   one of the character strings `"random"` (to select a seed uniformly at
#'   random) and `"minmax"` (to select a seed from the minmax set).
#' @param shuffle_data whether to first randomly shuffle the data.
#' @references De Silva, Vin, and Gunnar E. Carlsson. "Topological estimation
#'   using witness complexes." SPBG 4 (2004): 157-166.
#' @references Dłotko, Paweł. "Ball Mapper: A Shape Summary for Topological Data
#'   Analysis." (2019). Web.
NULL

#' @rdname landmarks_maxmin
#' @export
minmax_R <- function(
  x,
  dist_method = "euclidean"
) {

  # initialize minimum distance
  dist_min <- Inf
  # initialize minmax set
  mm_idx <- integer(0)

  # across all points
  for (idx in seq(nrow(x))) {

    # maximum distance from point
    dist_idx <- max(proxy::dist(x[idx, , drop = FALSE],
                                x,
                                method = dist_method))

    if (dist_idx == dist_min) {
      # if equal to reigning minimum distance, append to minmax set
      mm_idx <- c(mm_idx, idx)
    } else if (dist_idx < dist_min) {
      # if less than reigning minimum distance, replace and reinitialize
      dist_min <- dist_idx
      mm_idx <- c(idx)
    }

  }

  # return minmax subset
  mm_idx
}

#' @rdname landmarks_maxmin
#' @export
chebyshev_center_R <- minmax_R

#' @rdname landmarks_maxmin
#' @export
landmarks_maxmin <- function(
  x,
  dist_method = "euclidean",
  num_sets = NULL, radius = NULL,
  seed_index = 1L, shuffle_data = FALSE
) {
  stopifnot(is.matrix(x))
  stopifnot(seed_index >= 1L && seed_index <= nrow(x))
  # must specify a number of balls or a radius
  stopifnot(! is.null(num_sets) || ! is.null(radius))

  shuffle_idx <- NA
  if (shuffle_data){
    shuffle_idx <- sample(seq(nrow(x)))
    x <- x[, , drop = FALSE]
  }

  if (!is.null(num_sets)) {
    if (missing(dist_method) || toupper(dist_method) == "EUCLIDEAN") {
      lmk_idx <- landmark_maxmin(x, num_sets, seed_index - 1L)
    } else if (requireNamespace("proxy", quietly = TRUE)) {
      stopifnot(toupper(dist_method) %in%
                  toupper(proxy::pr_DB$get_entry_names()))
      lmk_idx <- vector(mode="integer", num_sets)
      lmk_idx[1L] <- seed_index
      lmk_min_dist <- rep(Inf, nrow(x))
      for (i in 2L:num_sets) {
        lmk_dist <- proxy::dist(x, x[lmk_idx[i - 1L], , drop = FALSE],
                                     method = dist_method)
        lmk_min_dist <- pmin(lmk_dist, lmk_min_dist)
        potential_idx <- setdiff(seq(nrow(x)), lmk_idx[c(1:i)])
        lmk_idx[i] <- potential_idx[which.max(lmk_min_dist[potential_idx])]
      }
    } else {
      stop(sprintf("Unsupported distance method passed: %s\n", dist_method))
    }
  } else if(! is.null(radius)) {
    stopifnot(toupper(dist_method) %in% toupper(proxy::pr_DB$get_entry_names()))
    f_dim <- ncol(x)
    f_size <- nrow(x)

    # algorithm and variable names as specified in Dlotko paper
    C = c() # create a list to store indices of centers/landmarks
    next_pt = seed_index # first landmark should be the seed point
    while(TRUE){
      C = append(C, next_pt) # add new point to list of landmarks

      # compute distance between landmark set and each point in the space
      dists = proxy::dist(matrix(x[C, ], ncol = f_dim), x, method = dist_method)
      sortedDists = matrix(apply(dists, 2, sort), ncol = f_size)

      # the next landmark is the point with greatest distance from the current
      # landmark set
      next_pt = which.max(sortedDists[1, ])
      # done when this max distance is < radius, i.e. when all pts are contained
      # in an radius-ball
      d = sortedDists[1, next_pt]
      if(d < radius) { break }
    }
    lmk_idx = C
  }
  if (is.na(shuffle_idx)) { lmk_idx } else { shuffle_idx[lmk_idx] }
}

#' @rdname landmarks_maxmin
#' @export
landmarks_maxmin_R <- function(
  x,
  dist_method = "euclidean", pick_method = "first",
  num_sets = NULL, radius = NULL, frac = FALSE,
  seed_index = 1L
) {
  # validate inputs
  stopifnot(is.matrix(x))
  stopifnot(toupper(dist_method) %in% toupper(proxy::pr_DB$get_entry_names()))

  # handle seed selection
  if (is.character(seed_index)) {
    seed_index <- switch (
      match.arg(seed_index, c("random", "minmax")),
      random = sample(nrow(x), size = 1L),
      minmax = sample(minmax_R(x, dist_method = dist_method), size = 1L)
    )
  }
  stopifnot(seed_index >= 1L && seed_index <= nrow(x))

  # initialize maxmin, free, and landmark index sets
  mm_idx <- seed_index
  lmk_idx <- vector(mode = "integer", nrow(x))
  free_idx <- seq(nrow(x))
  # strike seed and (other) duplicates index from free indices
  perm_idx <- c(mm_idx, free_idx[-mm_idx])
  free_idx[perm_idx[duplicated(x[perm_idx])]] <- 0L
  lmk_dist <- rep(Inf, times = nrow(x))

  # apply `frac` to `radius`
  if (frac) {
    # -+- a `chull()` function for arbitrary dimensions would expedite this -+-
    diameter <- max(proxy::dist(x))
    radius <- radius * diameter
  }

  for (i in seq_along(free_idx)) {

    # update vector of landmark points
    lmk_idx[[i]] <- mm_idx[[switch (
      match.arg(pick_method, c("first", "last", "random")),
      first = 1L,
      last = length(mm_idx),
      random = sample(length(mm_idx), 1L)
    )]]

    # update vector of free points
    if (free_idx[[lmk_idx[[i]]]] == 0L)
      stop("A duplicate landmark point was selected, in error.")
    free_idx[[lmk_idx[[i]]]] <- 0L

    # update landmark distances with distances from new landmark point
    lmk_dist <- pmin(lmk_dist,
                     proxy::dist(x[lmk_idx[[i]], , drop = FALSE],
                                 x,
                                 method = dist_method))
    # refresh the minimum radius
    min_rad <- max(lmk_dist)

    # exhaustion breaks
    if (all(free_idx == 0L)) break
    # parameter breaks
    if (! is.null(num_sets)) {
      # discontinue if the desired number of sets has been reached
      if (i >= num_sets) {
        # discontinue if there is no radius requirement
        if (is.null(radius)) {
          break
        } else {
          # continue if desired radius requires a greater number of sets
          if (min_rad <= radius) break
        }
      }
    } else {
      if (is.null(radius)) {
        # if neither parameter is specified, limit the set to 24 landmarks
        if (i >= 24L) break
      } else {
        # continue if desired radius requires more sets
        if (min_rad <= radius) break
      }
    }

    # obtain the maxmin subset
    mm_idx <- free_idx[free_idx != 0L]
    mm_idx <- mm_idx[lmk_dist[mm_idx] == min_rad]

  }

  # print warnings if a parameter was adjusted
  if (! is.null(num_sets)) {
    if (i > num_sets) {
      warning("Cover required ", i, " (> num_sets = ", num_sets, ") ",
              "balls of radius ", radius, ".")
    } else if (i < num_sets) {
      warning("Only ", i, " (< num_sets = ", num_sets, ") ",
              "distinct landmark points were found.")
    }
  }

  lmk_idx[seq(i)]
}
