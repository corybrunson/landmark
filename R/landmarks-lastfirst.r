#' @name landmarks_lastfirst
#' @title Neighborhood-based Landmark Sets
#' @author Jason Cory Brunson
#' @author Yara Skaf
#' @description Compute landmark sets based on nearest neighborhoods.
#' @details These functions adapt the maxmin procedure to produce landmark
#'   points dispersed according to the orders in which they are reached from
#'   each other, rather than to their distances from each other. (Say more.)
#'
#'   One, both, or neither of `num` and `cardinality` may be passed values. If
#'   neither is specified, then `num` is defaulted to the minimum of `24L` and
#'   the number of distinct rows of `x`. If the values yield neighborhoods that
#'   do not cover `x`, then, effectively, `num` is increased until the
#'   cardinality necessary to cover `x` is at most `cardinality`. To generte a
#'   complete landmark set, use `cardinality = 1L`.
#' @param x a data matrix.
#' @param y a data matrix of the same dimension as `x`; if `NULL`, taken to be
#'   `x`.
#' @param dist_method a character string specifying the distance metric to use;
#'   passed to `proxy::dist(method)`. Any distance measure in the \code{proxy}
#'   package is supported.
#' @param ties_method a character string specifying the method for handling
#'   ties; passed to `rank(ties.method)`. Only `"min"` and `"max"` have been
#'   tested and are recommended.
#' @param pick_method a character string specifying the method for selecting one
#'   among indistinguishable points, either `"first"` (the default), `"last"`,
#'   or `"random"`.
#' @param num a positive integer; the desired number of landmark points, or of
#'   sets in a neighborhood cover.
#' @param cardinality a positive integer; the desired cardinality of each
#'   landmark neighborhood, or of each set in a landmark cover.
#' @param frac logical; whether to treat `cardinality` as a fraction of the
#'   cardinality of `x`.
#' @param seed_index an integer (the first landmark to seed the algorithm) or
#'   one of the character strings `"random"` (to select a seed uniformly at
#'   random) and `"firstlast"` (to select a seed from the firstlast set).
#' @param engine character string specifying the implementation to use; one of
#'   `"C++"` or `"R"`. When not specified, the R engine is used.
#' @param cover logical; whether to include a column of cover sets (by member
#'   index) alongside columns of landmark indices and of minimum neighborhood
#'   cardinalities in a data frame.
#' @param tower logical; whether to include a column of covers (stored as
#'   maximal simplices of their nerves) alongside columns of landmark indices
#'   and of minimum neighborhood cardinalities in a data frame.
#' @param extend_num,extend_cardinality length-two numeric vectors used to
#'   extend landmark parameters for cover set construction. See [extension()].
#' @example inst/examples/ex-landmarks-lastfirst.r
NULL

#' @rdname landmarks_lastfirst
#' @export
firstlast <- function(
    x, y = NULL,
    dist_method = "euclidean", ties_method = "min"
) {
  
  # use distances from `x` if `y` is not specified
  if (is.null(y)) y <- x
  
  # update if/when C++ implementation is available
  firstlast_R(x = x, y = y,
              dist_method = dist_method, ties_method = ties_method)
}

firstlast_R <- function(
    x, y,
    dist_method = "euclidean", ties_method = "min"
) {
  
  # initialize colex-minimum rank-out-distance sequence
  seq_min <- rep(nrow(y), nrow(y))
  # initialize firstlast set
  fl_idx <- integer(0)
  
  # across all points
  for (idx in seq(nrow(x))) {
    
    # out-rank-distance sequence
    seq_idx <- sort(rank(proxy::dist(
      x[idx, , drop = FALSE],
      y,
      method = dist_method
    ), ties.method = ties_method))
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
lastfirst <- function(
    x, y = NULL,
    dist_method = "euclidean", ties_method = "min"
) {
  
  # use distances from `x` if `y` is not specified
  if (is.null(y)) {
    y <- x
    self <- TRUE
  } else {
    self <- FALSE
  }
  
  # update if/when C++ implementation is available
  lastfirst_R(x = x, y = y, self = self,
              dist_method = dist_method, ties_method = ties_method)
}

lastfirst_R <- function(
    x, y, self,
    dist_method = "euclidean", ties_method = "min"
) {
  
  # initialize revlex-maximum rank-in-distance sequence
  seq_max <- rep(nrow(y), nrow(y))
  # initialize lastfirst set
  lf_idx <- integer(0)
  
  # distance matrix
  dist_mat <- proxy::dist(
    x,
    y,
    method = dist_method
  )
  if (self) {
    # exclude diagonal
    dist_mat <- t(matrix(
      dist_mat[seq(nrow(x) * nrow(y)) %% (nrow(y) + 1L) != 1L],
      nrow = nrow(y) - 1L
    ))
  }
  # all in-rank-distance sequences
  rank_mat <- t(apply(
    t(apply(dist_mat, 1L, rank, ties.method = ties_method)),
    1L, sort
  ))
  # obtain the lastfirst subset
  lf_idx <- seq(nrow(rank_mat))
  for (j in seq(ncol(rank_mat))) {
    # points with maximum revlex rank-in-distance counts
    # = points with minimum lex rank-in-distance counts
    # = points with maximum lex rank-in-distance sequence
    lf_idx <- lf_idx[rank_mat[lf_idx, j] == max(rank_mat[lf_idx, j])]
    if (length(lf_idx) == 1L) break
  }
  
  # return firstlast subset
  lf_idx
}

#' @rdname landmarks_lastfirst
#' @export
landmarks_lastfirst <- function(
    x,
    dist_method = "euclidean", ties_method = "min", pick_method = "first",
    num = NULL, cardinality = NULL, frac = FALSE,
    seed_index = 1L,
    engine = NULL,
    cover = FALSE, tower = FALSE,
    extend_num = extension(mult = 0, add = 0),
    extend_cardinality = extension(mult = 0, add = 0)
) {
  # validate inputs
  stopifnot(is.matrix(x))
  dist_method <- tolower(dist_method)
  stopifnot(dist_method %in% tolower(proxy::pr_DB$get_entry_names()))
  pick_method <- match.arg(pick_method, c("first", "last", "random"))
  if (is.null(engine)) engine <- "R"
  engine <- match.arg(engine, c("C++", "R"))
  if (engine == "C++" && dist_method != "euclidean") {
    warning("C++ engine is available only for Euclidean distances; ",
            "using R engine instead.")
    engine <- "R"
  }
  if (engine == "R") {
    mult_num <- extend_num[[1L]]
    add_num <- extend_num[[2L]]
    mult_cardinality <- extend_cardinality[[1L]]
    add_cardinality <- extend_cardinality[[2L]]
  }
  
  # if neither parameter is specified, limit the set to 24 landmarks
  if (is.null(num) && is.null(cardinality)) {
    num <- min(nrow(unique(x)), 24L)
  }
  # apply `frac` to `cardinality`
  if (frac) {
    cardinality <- as.integer(max(1, cardinality * nrow(x)))
  }
  # validate parameters
  if (! is.null(num)) {
    num <- as.integer(num)
    if (is.na(num) || num < 1L || num > nrow(x))
      stop("`num` must be a positive integer and at most `nrow(x)`.")
  }
  if (! is.null(cardinality)) {
    cardinality <- as.integer(cardinality)
    if (is.na(cardinality) || cardinality < 1L || cardinality > nrow(x))
      stop("`cardinality` must be a positive integer and at most `nrow(x)`.")
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
      match.arg(seed_index, c("random", "firstlast")),
      random = sample(nrow(x), size = 1L),
      firstlast = {
        fl_idx <- firstlast(x,
                            dist_method = dist_method,
                            ties_method = ties_method)
        fl_idx[[1L]]
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
  stopifnot(seed_index >= 1L, seed_index <= nrow(x))
  
  # dispatch to implementations
  res <- switch (
    engine,
    `C++` = landmarks_lastfirst_cpp(
      x = x,
      num = if (is.null(num)) 0L else num,
      cardinality = if (is.null(cardinality)) 0L else cardinality,
      seed_index = seed_index, cover = cover
    ),
    R = landmarks_lastfirst_R(
      x = x,
      dist_method = dist_method,
      ties_method = ties_method,
      num = num, cardinality = cardinality,
      seed_index = seed_index, cover = cover, tower = tower,
      mult_num = mult_num, add_num = add_num,
      mult_cardinality = mult_cardinality,
      add_cardinality = add_cardinality
    )
  )
  
  # format list as a data frame
  # TODO: Be readier for edits to C++ implementations.
  stopifnot(is.list(res))
  if (length(res) == 1L) {
    res <- res[[1L]]
  } else {
    names(res) <- if (length(res) == 2L) {
      c("landmark", "cover_set")
    } else {
      c("landmark", "cardinality", "cover_set", "tower_max")
    }
    res <- res[! vapply(res, is.null, TRUE)]
    list_cols <- which(vapply(res, is.list, TRUE))
    for (list_col in list_cols) res[[list_col]] <- I(res[[list_col]])
    res <- as.data.frame(res)
  }
  
  # correct for permutation
  if (pick_method != "first") {
    if (is.list(res)) {
      res[["landmark"]] <- shuffle_idx[res[["landmark"]]]
      if (! is.null(res[["cover_set"]])) res[["cover_set"]] <- lapply(
        res[["cover_set"]],
        function(set) shuffle_idx[set]
      )
      if (! is.null(res[["tower_max"]])) {
        for (tower_i in seq_along(res[["tower_max"]])) {
          res[["tower_max"]][[tower_i]] <- lapply(
            res[["tower_max"]][[tower_i]],
            function(set) shuffle_idx[set]
          )
        }
      }
    } else {
      res <- shuffle_idx[res]
    }
  }
  
  # print warnings if a parameter was adjusted
  ext_num <- num * (1 + extend_num[[1L]]) + extend_num[[2L]]
  if (! is.null(num)) {
    if (NROW(res) > ext_num) {
      warning("Required ", NROW(res),
              " (> num = ", ext_num, ") ",
              "sets of cardinality ", cardinality, ".")
    } else if (NROW(res) < ext_num) {
      warning("Only ", NROW(res),
              " (< num = ", ext_num, ") ",
              "distinct landmark points were found.")
    }
  }
  
  # return landmarks
  res
}

landmarks_lastfirst_R <- function(
    x,
    dist_method = "euclidean", ties_method = "min",
    num = NULL, cardinality = NULL,
    seed_index = 1L, cover = FALSE, tower = FALSE,
    mult_num = 0, add_num = 0, mult_cardinality = 0, add_cardinality = 0
) {
  
  # initialize lastfirst, free, and landmark index sets
  lf_idx <- seed_index
  free_idx <- seq(nrow(x))
  lmk_idx <- vector(mode = "integer", nrow(x))
  # strike seed and (other) duplicates index from free indices
  # -+- assumes data have been permuted according to `pick_method` -+-
  perm_idx <- c(lf_idx, free_idx[-lf_idx])
  free_idx[perm_idx[duplicated(x[perm_idx, , drop = FALSE])]] <- 0L
  # initialize rank matrix and membership list
  lmk_rank <- matrix(NA, nrow = nrow(x), ncol = 0)
  # initialize minimum cardinality and associated number of sets to cover `x`
  cover_card <- nrow(x)
  cover_num <- 0L
  card_idx <- vector(mode = "integer", nrow(x))
  if (cover || tower) cover_idx <- list()
  if (tower) tower_idx <- list()
  
  # recursively construct landmark set
  for (i in seq_along(free_idx)) {
    
    # update vector of landmark points
    # -+- assumes data have been permuted according to `pick_method` -+-
    lmk_idx[[i]] <- lf_idx[[1L]]
    
    # update vector of free points
    if (free_idx[[lmk_idx[[i]]]] == 0L)
      stop("A duplicate landmark point was selected, in error.")
    free_idx[[lmk_idx[[i]]]] <- 0L
    
    # augment in-ranks from new landmark point
    lmk_rank <- cbind(
      # each row contains the in-ranks from previous landmark points
      lmk_rank,
      # in-ranks of all points from newest landmark point
      rank(proxy::dist(
        x[lmk_idx[[i]], , drop = FALSE],
        x,
        method = dist_method), ties.method = ties_method)
    )
    
    # minimum in-rank from landmarks to `x`
    min_rank <- max(pmin(lmk_rank[c(free_idx, lmk_idx), 1L],
                         lmk_rank[c(free_idx, lmk_idx), ncol(lmk_rank)]))
    # update the minimum cardinality necessary to cover `x`
    if (is.null(cardinality) && ! is.null(num) && i <= num)
      cover_card <- min_rank
    # update radius vector
    card_idx[[i]] <- cover_card
    # update membership list
    if (cover || tower) {
      cover_idx <- if (is.null(cardinality)) {
        # -+- will need to parse later -+-
        wh_idx <- which(lmk_rank[, ncol(lmk_rank)] <=
                          cover_card * (1 + mult_cardinality) + add_cardinality)
        c(cover_idx,
          list(cbind(idx = wh_idx, rank = lmk_rank[wh_idx, ncol(lmk_rank)])))
      } else {
        # -+- will not need to parse later -+-
        c(cover_idx, list(which(lmk_rank[, ncol(lmk_rank)] <=
                                  cardinality * (1 + mult_cardinality) +
                                  add_cardinality)))
      }
    }
    # update nerve list (stored as lists of maximal simplices)
    if (tower) {
      cover_sets <- lapply(cover_idx, function(mat) {
        unname(mat[mat[, "rank"] <=
                     cover_card * (1 + mult_cardinality) +
                     add_cardinality, "idx"])
      })
      cover_mems <- transpose_list(cover_sets)
      nerve <-
        simplextree::nerve(simplextree::simplex_tree(), cover_mems, k = 2L)
      nerve <-
        simplextree::ltraverse(simplextree::maximal(nerve), f = as.integer)
      tower_idx <- c(tower_idx, list(nerve))
    }
    
    # exhaustion breaks
    if (all(free_idx == 0L)) break
    # update the minimum number of sets necessary to cover `x`
    if (is.null(num) && ! is.null(cardinality) && min_rank > cardinality)
      cover_num <- i + 1L
    # parameter breaks
    if ((is.null(num) || i >= num * (1L + mult_num) + add_num) &&
        (cover_num == 0L || i >= cover_num * (1L + mult_num) + add_num) &&
        (is.null(cardinality) || min_rank <= cardinality)) break
    
    # sort each available point's in-ranks to the landmark points
    lmk_rank[] <- t(apply(lmk_rank, 1L, sort))
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
  
  # return data
  list(
    # restrict to selected landmarks
    lmk_idx[seq(i)], card_idx[seq(i)],
    if (cover) {
      # parse extraneous members
      lapply(cover_idx, function(mat) {
        unname(mat[mat[, "rank"] <=
                     cover_card * (1 + mult_cardinality) +
                     add_cardinality, "idx"])
      })
    },
    if (tower) tower_idx
  )
}
