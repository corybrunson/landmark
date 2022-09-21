#' @name transpose_list
#' @title Convert between lists of members and lists of memberships
#' @description Covers may be stored as lists of set members or as lists of
#'   element memberships. This function converts one to the other (the operation
#'   is an involution).
#' @param x A list of positive integer vectors. Names will be ignored.
#' @param sparse Logical; whether to use a slower `lapply()` implementation
#'   rather than a speedier but possibly memory-intensive implementation that
#'   constructs the element-set matrix. (Neither implementation currentnly takes
#'   advantage of sparse matrix objects.)
#' @example inst/examples/ex-transpose-list.r
#' @export
transpose_list <- function(x, sparse = TRUE) {
  # ensure that `x` is a list of positive integer vectors
  stopifnot(
    is.list(x),
    all(vapply(x, is.vector, FALSE)),
    all(vapply(x, is.numeric, FALSE)),
    all(vapply(x, function(vec) all(vec > 0 & vec %% 1 == 0), FALSE))
  )
  # vectors of elements and of sets (or vice-versa)
  elts <- unique(unlist(x))
  elts <- seq(1L, max(elts))
  sets <- seq_along(x)
  # 'transpose' list of element vectors to list of set vectors
  if (sparse) {
    lapply(
      elts,
      function(elt) sets[which(vapply(x, function(vec) elt %in% vec, FALSE))]
    )
  } else {
    elt_by_set <- sapply(x, function(vec) elts %in% vec)
    apply(elt_by_set, 1L, function(row) sets[row])
  }
}
