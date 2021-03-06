% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/landmarks-maxmin.r
\name{landmarks_maxmin}
\alias{landmarks_maxmin}
\alias{minmax}
\alias{maxmin}
\title{Ball-based Landmark Sets}
\usage{
minmax(x, y = NULL, dist_method = "euclidean")

maxmin(x, y = NULL, dist_method = "euclidean")

landmarks_maxmin(
  x,
  dist_method = "euclidean",
  pick_method = "first",
  num = NULL,
  radius = NULL,
  frac = FALSE,
  seed_index = 1L,
  engine = NULL,
  cover = FALSE,
  extend_num = extension(mult = 0, add = 0),
  extend_radius = extension(mult = 0, add = 0)
)
}
\arguments{
\item{x}{a data matrix.}

\item{y}{a data matrix of the same dimension as \code{x}; if \code{NULL}, taken to be
\code{x}.}

\item{dist_method}{a character string specifying the distance metric to use;
passed to \code{proxy::dist(method)}. Any distance measure in the \code{proxy}
package is supported.}

\item{pick_method}{a character string specifying the method for selecting one
among indistinguishable points, either \code{"first"} (the default), \code{"last"},
or \code{"random"}.}

\item{num}{a positive integer; the desired number of landmark points, or of
sets in a ball cover.}

\item{radius}{a positive number; the desired radius of each landmark ball, or
of each set in a ball cover.}

\item{frac}{logical; whether to treat \code{radius} as a fraction of the diameter
of \code{x}.}

\item{seed_index}{an integer (the first landmark to seed the algorithm) or
one of the character strings \code{"random"} (to select a seed uniformly at
random) and \code{"minmax"} (to select a seed from the minmax set).}

\item{engine}{character string specifying the implementation to use; one of
\code{"original"}, \code{"C++"}, or \code{"R"}. When not specified, the R engine is used.}

\item{cover}{logical; whether to return a data frame of landmark indices and
cover sets (by member index) rather than only a vector of landmark indices.}

\item{extend_num, extend_radius}{length-two numeric vectors used to extend
landmark parameters for cover set construction. See \code{\link[=extension]{extension()}}.}
}
\description{
Compute landmark sets based on fixed-radius balls.
}
\details{
This function uses the maxmin procedure to produce a set of evenly
spaced landmark points from a data set. Maxmin is a simple greedy algorithm
that is relatively efficient, but it has a tendency to pick out extremal
points.

One, both, or neither of \code{num} and \code{radius} may be passed values. If
neither is specified, then \code{num} is defaulted to \code{24L}. If the values yield
balls that do not cover \code{x}, then their number is increased until the
radius necessary to cover \code{x} is at most \code{radius}. To generte a complete
landmark set, use \code{radius = 0L}.
}
\examples{
set.seed(3)
# small circle sample
X <- tdaunif::sample_circle(n = 12L)
# random seed index
l <- landmarks_maxmin(X, seed_index = "random")
# plot landmark order at point positions
plot(X, asp = 1, pch = NA)
text(X, labels = order(l))
# minmax seed index
l <- landmarks_maxmin(X, seed_index = "minmax")
# plot landmark order at point positions
plot(X, asp = 1, pch = NA)
text(X, labels = order(l))
}
\references{
De Silva, Vin, and Gunnar E. Carlsson. "Topological estimation
using witness complexes." SPBG 4 (2004): 157-166.

Dłotko, Paweł. "Ball Mapper: A Shape Summary for Topological Data
Analysis." (2019). Web.
}
\author{
Matt Piekenbrock

Jason Cory Brunson

Yara Skaf
}
