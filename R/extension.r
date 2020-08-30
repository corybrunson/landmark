#' @name extension
#' @title Generate and validate an extension vector for landmark parameters
#' @description This is a convenience function, modeled after
#'   [ggplot2::expansion()], for generating extension vectors for the `extend_*`
#'   parameters of [landmarks_maxmin()] and [landmarks_lastfirst()].
#' @details The parameter `num` is extended to `num * (1 + extend_num[[1L]]) +
#'   extend_num[[2L]]`, and other parameters are extended analogously.
#'
#'   These extensions are used to generate covers whose sets are more numerous
#'   or larger than strictly necessary to cover the finite metric space. When
#'   performing topological analysis, e.g. constructing a nerve complex,
#'   redundant overlap may be necessary to capture essential features of the
#'   data.
#' @param mult multiplicative extension factor; added to `1` and multiplied by
#'   the parameter.
#' @param add additive extension factor; added to parameter after
#'   multiplication.
#' @export
extension <- function(mult = 0, add = 0) {
  if (! (is.numeric(mult) && (length(mult) == 1L) && mult >= 0 &&
         is.numeric(add) && (length(add) == 1L) && add >= 0)) {
    stop("`mult` and `add` must be non-negative scalars.")
  }
  c(mult, add)
}
