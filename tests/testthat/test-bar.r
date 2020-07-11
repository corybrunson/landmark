context("interval with different densities close to endpoints")

peg <- matrix(c(1, 2, 4, 4))
bar <- matrix(c(-1, -.5, 0, .75, .875, 1))

# maxmin landmarks in C++

test_that("balls contain points exactly radius from center", {
  mm_r1 <- landmarks_maxmin(bar, radius = 1, seed_index = 3L)
  expect_identical(mm_r1, 3L)
})

test_that("2-set cover uses endpoints", {
  mm_n2 <- landmarks_maxmin(bar, num_sets = 2L, engine = "C++")
  expect_identical(mm_n2, c(1L, 6L))
})

test_that("half-diameter-radius cover uses endpoints only", {
  mm_r1 <- landmarks_maxmin(bar, radius = 1, engine = "C++")
  expect_identical(mm_r1, c(1L, 6L))
})

test_that("complete landmark set grows leftward before righward", {
  mm_all <- landmarks_maxmin(bar, num_sets = 6L, engine = "C++")
  expect_identical(mm_all, c(1L, 6L, 3L, 2L, 4L, 5L))
})

# maxmin landmarks in R

test_that("2-set cover uses endpoints", {
  mm_n2 <- landmarks_maxmin(bar, num_sets = 2L, engine = "R")
  expect_identical(mm_n2, c(1L, 6L))
})

test_that("half-diameter-radius cover uses endpoints only", {
  mm_r1 <- landmarks_maxmin(bar, radius = 1, engine = "R")
  expect_identical(mm_r1, c(1L, 6L))
})

test_that("complete landmark set grows leftward before righward", {
  mm_all <- landmarks_maxmin(bar, num_sets = 6L, engine = "R")
  expect_identical(mm_all, c(1L, 6L, 3L, 2L, 4L, 5L))
})

# firstlast chebycenter

test_that("firstlast procedure prefers marginalizes multiple-points", {
  # both singleton points are centers because they are not the double-point
  expect_identical(firstlast(peg), c(1L, 2L))
  # only zero is a center because it is equidistant from the endpoints
  expect_identical(firstlast(bar), c(3L))
})

# lastfirst landmarks

test_that("landmark set is obtained starting from duplicate point", {
  expect_error(landmarks_lastfirst(peg, seed_index = 3L), NA)
})

test_that("invalid numbers of sets prompt warnings", {
  expect_warning(landmarks_lastfirst(peg, num_sets = 1L, cardinality = 2L),
                 "cardinality")
  expect_warning(landmarks_lastfirst(peg, num_sets = 4L),
                 "landmark")
})

test_that("number of cardinality-3 cover sets of peg depends on seed", {
  lf_k3_s1 <- landmarks_lastfirst(peg, cardinality = 3L, seed_index = 1L)
  expect_identical(lf_k3_s1, c(1L))
  lf_k3_s2 <- landmarks_lastfirst(peg, cardinality = 3L, seed_index = 2L)
  expect_identical(lf_k3_s2, c(2L))
  lf_k3_s3 <- landmarks_lastfirst(peg, cardinality = 3L, seed_index = 3L)
  expect_identical(lf_k3_s3, c(3L, 1L))
})

test_that("cardinality-2 cover of peg requires two points", {
  lf_k2_s1 <- landmarks_lastfirst(peg, cardinality = 2L, seed_index = 1L)
  expect_identical(lf_k2_s1, c(1L, 3L))
  lf_k2_s2 <- landmarks_lastfirst(peg, cardinality = 2L, seed_index = 2L)
  expect_identical(lf_k2_s2, c(2L, 3L))
  lf_k2_s3 <- landmarks_lastfirst(peg, cardinality = 2L, seed_index = 3L)
  expect_identical(lf_k2_s3, c(3L, 1L))
})

test_that("complete landmark set of peg grows left-to-right", {
  lf_all <- suppressWarnings(landmarks_lastfirst(peg, num_sets = 4L,
                                                 seed_index = 1L))
  expect_identical(lf_all, c(1L, 3L, 2L))
})

test_that("2-set cover of bar uses endpoints", {
  lf_n2 <- landmarks_lastfirst(bar, num_sets = 2L)
  expect_identical(lf_n2, c(1L, 6L))
})

test_that("cardinality-3 cover of bar uses endpoints only", {
  lf_k3 <- landmarks_lastfirst(bar, cardinality = 3L)
  expect_identical(lf_k3, c(1L, 6L))
})

test_that("cardinality-2 cover of bar uses end- and median points", {
  lf_k2 <- landmarks_lastfirst(bar, cardinality = 2L)
  expect_identical(lf_k2, c(1L, 6L, 3L, 4L))
})

test_that("complete landmark set of bar grows left-to-right", {
  lf_all <- landmarks_lastfirst(bar, num_sets = 6L)
  expect_identical(lf_all, c(1L, 6L, 3L, 4L, 2L, 5L))
})
