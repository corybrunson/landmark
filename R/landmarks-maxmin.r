#' @name landmarks_maxmin
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
#' @references De Silva, Vin, and Gunnar E. Carlsson. "Topological estimation
#'   using witness complexes." SPBG 4 (2004): 157-166.
#' @references Dłotko, Paweł. "Ball Mapper: A Shape Summary for Topological Data
#'   Analysis." (2019). Web.
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
#' @param cover logical; whether to return a data frame of landmark indices and
#'   cover sets (by member index) rather than only a vector of landmark indices.
#' @example inst/examples/ex-landmarks-maxmin.r
NULL

#' @rdname landmarks_maxmin
#' @export
minmax <- function(
  x,
  dist_method = "euclidean"
) {
  # update if/when C++ implementation is available
  minmax_R(x = x, dist_method = dist_method)
}

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
  dist_method = "euclidean", pick_method = "first",
  num_sets = NULL, radius = NULL, frac = FALSE,
  seed_index = 1L, cover = FALSE,
  engine = NULL
) {
  # validate inputs
  stopifnot(is.matrix(x))
  dist_method <- tolower(dist_method)
  stopifnot(dist_method %in% tolower(proxy::pr_DB$get_entry_names()))
  pick_method <- match.arg(pick_method, c("first", "last", "random"))
  if (is.null(engine)) engine <- "R"
  engine <- match.arg(engine, c("original", "C++", "R"))
  if (engine == "C++" && dist_method != "euclidean")
    warning("C++ engine is available only for Euclidean distances; ",
            "using R engine instead.")

  # if neither parameter is specified, limit the set to 24 landmarks
  if (is.null(num_sets) && is.null(radius)) {
    num_sets <- min(nrow(unique(x)), 24L)
  }
  # apply `frac` to `radius`
  if (frac) {
    radius <- max(1, radius * nrow(x))
  }
  # validate parameters
  if (! is.null(num_sets)) {
    num_sets <- as.integer(num_sets)
    if (is.na(num_sets) || num_sets < 1L || num_sets > nrow(x))
      stop("`num_sets` must be a positive integer and at most `nrow(x)`.")
  }
  if (! is.null(radius)) {
    if (is.na(radius) || radius <= 0 || radius == Inf)
      stop("`radius` must be a finite non-negative number.")
  }

  # permute rows of `x` according to `pick_method`
  if (pick_method != "first") {
    shuffle_idx <- switch (
      pick_method,
      first = seq(nrow(x)),
      last = seq(nrow(x), 1L),
      random = sample(nrow(x))
    )
    x <- x[shuffle_idx, , drop = FALSE]
  }
  # handle seed selection
  if (is.character(seed_index)) {
    seed_index <- switch (
      match.arg(seed_index, c("random", "minmax")),
      random = sample(nrow(x), size = 1L),
      minmax = {
        mm_idx <- minmax(x,
                         dist_method = dist_method)
        mm_idx[[1L]]
      }
    )
  } else {
    # reset input seed index accordingly
    if (pick_method != "first") {
      seed_index <- switch (
        pick_method,
        first = seed_index,
        last = nrow(x) + 1L - seed_index,
        random = which(shuffle_idx == seed_index)
      )
    }
  }
  stopifnot(seed_index >= 1L && seed_index <= nrow(x))

  # dispatch to implementations
  res <- switch (
    engine,
    original = landmarks_maxmin_orig(
      x = x,
      dist_method = dist_method,
      num_sets = num_sets, radius = radius,
      seed_index = seed_index
    ),
    `C++` = landmarks_maxmin_cpp(
      x = x,
      num_sets = if (is.null(num_sets)) 0L else num_sets,
      radius = if (is.null(radius)) -1L else radius,
      seed_index = seed_index, cover = cover
    ),
    R = landmarks_maxmin_R(
      x = x,
      dist_method = dist_method,
      num_sets = num_sets, radius = radius,
      seed_index = seed_index, cover = cover
    )
  )

  # format list as a data frame
  stopifnot(is.list(res))
  if (length(res) == 1L) {
    res <- res[[1]]
  } else {
    res <- data.frame(landmark = res[[1]], cover_set = I(res[[2]]))
  }
  # correct for permutation
  if (pick_method != "first") {
    if (is.list(res)) {
      res[[1]] <- shuffle_idx[res[[1]]]
      res[[2]] <- lapply(res[[2]], function(set) shuffle_idx[set])
    } else {
      res <- shuffle_idx[res]
    }
  }

  # print warnings if a parameter was adjusted
  if (! is.null(num_sets)) {
    if (NROW(res) > num_sets) {
      warning("Required ", NROW(res),
              " (> num_sets = ", num_sets, ") ",
              "sets of radius ", radius, ".")
    } else if (NROW(res) < num_sets) {
      warning("Only ", NROW(res),
              " (< num_sets = ", num_sets, ") ",
              "distinct landmark points were found.")
    }
  }

  # return landmarks
  res
}

#' @rdname landmarks_maxmin
#' @export
landmarks_maxmin_orig <- function(
  x,
  dist_method = "euclidean", pick_method = "first",
  num_sets = NULL, radius = NULL,
  seed_index = 1L
) {
  stopifnot(is.matrix(x))
  stopifnot(seed_index >= 1L && seed_index <= nrow(x))
  # must specify a number of balls or a radius
  stopifnot(! is.null(num_sets) || ! is.null(radius))

  shuffle_idx <- switch (
    pick_method,
    first = NA,
    last = rev(seq(nrow(x))),
    random = sample(seq(nrow(x)))
  )

  if (!is.null(num_sets)) {
    if (missing(dist_method) || toupper(dist_method) == "EUCLIDEAN") {
      lmk_idx <- landmark_maxmin(x, num_sets, seed_index)
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
  list(if (is.na(shuffle_idx)) { lmk_idx } else { shuffle_idx[lmk_idx] })
}

#' @rdname landmarks_maxmin
#' @export
landmarks_maxmin_R <- function(
  x,
  dist_method = "euclidean", pick_method = "first",
  num_sets = NULL, radius = NULL, frac = FALSE,
  seed_index = 1L, cover = FALSE
) {

  # initialize maxmin, free, and landmark index sets
  mm_idx <- seed_index
  lmk_idx <- vector(mode = "integer", nrow(x))
  free_idx <- seq(nrow(x))
  # strike seed and (other) duplicates index from free indices
  perm_idx <- c(mm_idx, free_idx[-mm_idx])
  free_idx[perm_idx[duplicated(x[perm_idx])]] <- 0L
  # initialize distance vector and membership list
  lmk_dist <- rep(Inf, times = nrow(x))
  min_rad <- Inf
  if (cover) cover_idx <- list()

  for (i in seq_along(free_idx)) {

    # update vector of landmark points
    # -+- assumes data have been permuted according to `pick_method` -+-
    lmk_idx[[i]] <- mm_idx[[1L]]

    # update vector of free points
    if (free_idx[[lmk_idx[[i]]]] == 0L)
      stop("A duplicate landmark point was selected, in error.")
    free_idx[[lmk_idx[[i]]]] <- 0L

    # update landmark distances with distances from new landmark point
    lmk_dist <- cbind(
      # minimum distances from previous landmark points
      lmk_dist,
      # distances of all points from newest landmark point
      proxy::dist(x[lmk_idx[[i]], , drop = FALSE],
                  x,
                  method = dist_method)[1, ])

    # refresh the minimum radius necessary to cover `x`
    min_rad <- max(pmin(lmk_dist[, 1L], lmk_dist[, 2L]))
    # update membership list
    if (cover) {
      cover_idx <- if (is.null(radius)) {
        # -+- will need to parse later -+-
        wh_idx <- which(lmk_dist[, 2L] <= min_rad)
        c(cover_idx,
          list(cbind(idx = wh_idx, dist = lmk_dist[wh_idx, 2L])))
      } else {
        # -+- will not need to parse later -+-
        c(cover_idx, list(which(lmk_dist[, 2L] <= radius)))
      }
    }

    # exhaustion breaks
    if (all(free_idx == 0L)) break
    # parameter breaks
    if ((is.null(num_sets) || i >= num_sets) &&
        (is.null(radius) || min_rad <= radius)) break

    # collapse distances to the minimum to each point
    lmk_dist <- pmin(lmk_dist[, 1L], lmk_dist[, 2L])
    # obtain the maxmin subset
    mm_idx <- free_idx[free_idx != 0L]
    mm_idx <- mm_idx[lmk_dist[mm_idx] == min_rad]

  }

  # restrict to selected landmarks
  lmk_idx <- lmk_idx[seq(i)]
  # return data
  if (cover) {
    if (is.null(radius)) {
      # parse extraneous members
      cover_idx <- lapply(cover_idx, function(mat) {
        unname(mat[mat[, "dist"] <= min_rad, "idx"])
      })
    }
    # return list of landmark indices and cover membership vectors
    list(lmk_idx, cover_idx)
  } else {
    # return vector of landmark indices
    list(lmk_idx)
  }
}
