#' @title Firstlast-generated Landmarks
#' @description Computes a landmark set based on nearest neighborhoods.
#' @details This function adapts the maxmin procedure to produce landmark points
#'   dispersed according to the orders in which they are reached from each
#'   other, rather than to their distances from each other. (Say more.)
#' @name landmarks_lastfirst
#' @param x a data matrix.
#' @param n the number of landmarks requested.
#' @param k the desired cardinality of each cover set.
#' @param dist_method the distance metric to use. Any distance measure in the
#'   \code{proxy} package is supported.
#' @param seed_index the first landmark to seed the algorithm.
#' @param shuffle_data whether to first randomly shuffle the data.
#' @export
landmarks_lastfirst_yara <- function(
  x, n=NULL, k=NULL,
  dist_method = "euclidean", seed_index = 1, shuffle_data=FALSE
){
  stopifnot(is.matrix(x))
  stopifnot(seed_index >= 1L && seed_index <= nrow(x))
  stopifnot(!is.null(n) || !is.null(eps)) # must specify a number of balls or a radius

  shuffle_idx <- NA
  if (shuffle_data){
    shuffle_idx <- sample(seq(nrow(x)))
    x <- x[,,drop=FALSE]
  }

  if(!is.null(n)){
    if (missing(dist_method) || toupper(dist_method) == "EUCLIDEAN"){
      landmark_idx <- landmark_maxmin(x, n, seed_index-1L)
    } else if (requireNamespace("proxy", quietly = TRUE)){
      stopifnot(toupper(dist_method) %in% toupper(proxy::pr_DB$get_entry_names()))
      landmark_idx <- vector(mode="integer", n)
      landmark_idx[1L] <- seed_index
      landmark_min_dist <- rep(Inf, nrow(x))
      for (i in 2L:n){
        landmark_dist <- proxy::dist(x, x[landmark_idx[i-1L],,drop=FALSE], method = dist_method)
        landmark_min_dist <- pmin(landmark_dist, landmark_min_dist)
        potential_idx <- setdiff(seq(nrow(x)), landmark_idx[c(1:i)])
        landmark_idx[i] <- potential_idx[which.max(landmark_min_dist[potential_idx])]
      }
    } else {
      stop(sprintf("Unsupported distance method passed: %s\n", dist_method))
    }
  } else if(!is.null(eps)){
    stopifnot(toupper(dist_method) %in% toupper(proxy::pr_DB$get_entry_names()))
    f_dim <- ncol(x)
    f_size <- nrow(x)

    # algorithm and variable names as specified in Dlotko paper
    C = c() # create a list to store indices of centers/landmarks
    next_pt = seed_index # first landmark should be the seed point
    while(TRUE){
      C = append(C, next_pt) # add new point to list of landmarks

      # compute distance between landmark set and each point in the space
      dists = proxy::dist(matrix(x[C,], ncol=f_dim), x, method = dist_method)
      sortedDists = matrix(apply(dists,2,sort),ncol=f_size)

      # the next landmark is the point with greatest distance from the current landmark set
      next_pt = which.max(sortedDists[1,])
      d = sortedDists[1,next_pt] # done when this max distance is < eps, i.e. when all pts are contained in an eps-ball
      if(d < eps){break}
    }
    landmark_idx = C
  }
  if (is.na(shuffle_idx)){ landmark_idx } else { shuffle_idx[landmark_idx] }
}
