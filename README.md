
<!-- README.md is generated from README.Rmd. Please edit that file -->

# landmark

Calculate landmark sets for finite metric spaces using the maxmin
procedure (for fixed-radius balls) or an adaptation of it for rank data
(for roughly fixed-cardinality nearest neighborhoods).

``` r
(x <- matrix(c(1, 2, 4, 4), dimnames = list(letters[1:4], "x")))
#>   x
#> a 1
#> b 2
#> c 4
#> d 4
```

## `maxmin` procedure

The original `maxmin` procedure produces a landmark set for covering a
point cloud with either of two minimal ball covers:

  - a minimum number of balls of fixed common radius
  - a fixed number of balls of minimum common radius

<!-- end list -->

``` r
landmarks_maxmin(x, eps = 3.5)
#> [1] 1
landmarks_maxmin(x, eps = 1.5)
#> [1] 1 3
landmarks_maxmin(x, eps = 0.5)
#> [1] 1 3 2
```

## “`lastfirst`” procedure

An adaptation of `maxmin` to ranked distances will produce a landmark
set for covering a point cloud with either of two minimal neighborhood
covers:

  - a minimum number of neighborhoods of roughly fixed common
    cardinality
  - a fixed number of neighborhoods of minimal roughly common
    cardinality

<!-- end list -->

``` r
landmarks_lastfirst_cory(x, cardinality = 4L, seed_index = 1L)
#> [1] 1
landmarks_lastfirst_cory(x, cardinality = 3L, seed_index = 1L)
#> [1] 1
landmarks_lastfirst_cory(x, cardinality = 2L, seed_index = 1L)
#> [1] 1 3
landmarks_lastfirst_cory(x, cardinality = 1L, seed_index = 1L)
#> [1] 1 3 2
```

# workflow

Foundations have been adapted from `peekxc/Mapper`.

We should each edit only our respective functions
`landmarks_lastfirst_*()` as we attempt to implement the algorithm
discussed on 2020-05-21, in the branches `lastfirst/*` (`cory` and
`yara`: experimental; `main`: consolidated).

If changes are made to other files then they should be committed to a
new branch first.

# reference

Rigorous definitions, algorithms, examples, and exposition are underway
at [this Overleaf project](https://www.overleaf.com/read/fpjrtgfjstyx).
