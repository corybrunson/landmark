
<!-- README.md is generated from README.Rmd. Please edit that file -->

# landmark

Calculate landmark sets for finite metric spaces using the maxmin
procedure (for fixed-radius balls) or an adaptation of it for rank data
(for roughly fixed-cardinality nearest
neighborhoods).

``` r
(x <- matrix(c(-1, -.5, 0, .75, .875, 1), dimnames = list(letters[1:6], "x")))
#>        x
#> a -1.000
#> b -0.500
#> c  0.000
#> d  0.750
#> e  0.875
#> f  1.000
plot(cbind(x, 0), asp = 1, pch = 16)
text(cbind(x, .05), labels = rownames(x))
```

<img src="man/figures/README-example-1.png" width="100%" />

## `maxmin` procedure

The original `maxmin` procedure produces a landmark set for covering a
point cloud with either of two minimal ball covers:

  - a minimum number of balls of fixed common radius
  - a fixed number of balls of minimum common radius

<!-- end list -->

``` r
x[landmarks_maxmin(x, eps = 1.5), , drop = FALSE]
#>    x
#> a -1
#> f  1
x[landmarks_maxmin(x, eps = 0.5), , drop = FALSE]
#>      x
#> a -1.0
#> f  1.0
#> c  0.0
#> b -0.5
x[landmarks_maxmin(x, eps = 0.25), , drop = FALSE]
#>       x
#> a -1.00
#> f  1.00
#> c  0.00
#> b -0.50
#> d  0.75
x[landmarks_maxmin(x, eps = 0.125), , drop = FALSE]
#>        x
#> a -1.000
#> f  1.000
#> c  0.000
#> b -0.500
#> d  0.750
#> e  0.875
```

## `lastfirst` procedure

An adaptation of `maxmin` to ranked distances will produce a landmark
set for covering a point cloud with either of two minimal neighborhood
covers:

  - a minimum number of neighborhoods of fixed common approximate
    cardinality
  - a fixed number of neighborhoods of minimal approximate cardinality

(Cardinality is only exact up to ties, which may be handled different
ways.)

``` r
x[landmarks_lastfirst_R(x, cardinality = 4L, seed_index = 1L), , drop = FALSE]
#>    x
#> a -1
#> f  1
x[landmarks_lastfirst_R(x, cardinality = 3L, seed_index = 1L), , drop = FALSE]
#>    x
#> a -1
#> f  1
x[landmarks_lastfirst_R(x, cardinality = 2L, seed_index = 1L), , drop = FALSE]
#>       x
#> a -1.00
#> f  1.00
#> c  0.00
#> d  0.75
x[landmarks_lastfirst_R(x, cardinality = 1L, seed_index = 1L), , drop = FALSE]
#>        x
#> a -1.000
#> f  1.000
#> c  0.000
#> d  0.750
#> b -0.500
#> e  0.875
```

# references

This package was spun off from [the Mapper
package](https://github.com/peekxc/Mapper/).

A rigorous mathematical treatment is underway at [this Overleaf
project](https://www.overleaf.com/read/fpjrtgfjstyx).
