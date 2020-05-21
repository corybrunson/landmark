
<!-- README.md is generated from README.Rmd. Please edit that file -->

# maxmin

Calculate landmark sets for finite metric spaces using the maxmin
procedure (for fixed-radius balls) or an adaptation of it for rank data
(for roughly fixed-cardinality nearest neighborhoods).

## `maxmin` procedure

The original `maxmin` procedure produces a landmark set for covering a
point cloud with either of two minimal ball covers:

  - a minimum number of balls of fixed common radius
  - a fixed number of balls of minimum common radius

## “`lastfirst`” procedure

An adaptation of `maxmin` to ranked distances will produce a landmark
set for covering a point cloud with either of two minimal neighborhood
covers:

  - a minimum number of neighborhoods of roughly fixed common
    cardinality
  - a fixed number of neighborhoods of minimal roughly common
    cardinality

# workflow

Foundations have been adapted from `peekxc/Mapper`, with `landmarks()`
renamed to `get_landmarks_maxmin()` to avoid variable overloading. We
should each edit only our respective functions
`get_landmarks_lastfirst_*()` as we attempt to implement the algorithm
discussed on 2020-05-21. If changes must be made to other files then
they should be committed to a new branch first.
