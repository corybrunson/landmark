context("interval with different densities close to endpoints")

peg <- matrix(c(1, 2, 4, 4))
bar <- matrix(c(-1, -.5, 0, .75, .875, 1))

# maxmin

test_that("minmax and maxmin procedures work on a single data matrix", {
  # both singleton points are centers because they are not the double-point
  expect_identical(minmax(bar), c(3L))
  # only zero is a center because it is equidistant from the endpoints
  expect_identical(maxmin(bar), c(1L, 2L, 3L))
})

test_that("minmax and maxmin procedures work on a multiple data matrices", {
  # nearest point(s) in `peg` to farthest point(s) in `bar`
  expect_identical(minmax(peg, bar), c(1L))
  # farthest point(s) in `peg` to nearest point(s) in `bar`
  expect_identical(maxmin(peg, bar), c(3L, 4L))
})

# maxmin landmarks in C++

test_that("balls contain points exactly radius from center", {
  mm_r1 <- landmarks_maxmin(bar, radius = 1, seed_index = 3L)
  expect_identical(mm_r1, 3L)
})

test_that("2-set cover uses endpoints", {
  mm_n2 <- landmarks_maxmin(bar, num = 2L, engine = "C++")
  expect_identical(mm_n2, c(1L, 6L))
})

test_that("half-diameter-radius cover uses endpoints only", {
  mm_r1 <- landmarks_maxmin(bar, radius = 1, engine = "C++")
  expect_identical(mm_r1, c(1L, 6L))
})

test_that("complete landmark set grows leftward before righward", {
  mm_all <- landmarks_maxmin(bar, num = 6L, engine = "C++")
  expect_identical(mm_all, c(1L, 6L, 3L, 2L, 4L, 5L))
})

# maxmin landmarks in R

test_that("2-set cover uses endpoints", {
  mm_n2 <- landmarks_maxmin(bar, num = 2L, engine = "R")
  expect_identical(mm_n2, c(1L, 6L))
})

test_that("half-diameter-radius cover uses endpoints only", {
  mm_r1 <- landmarks_maxmin(bar, radius = 1, engine = "R")
  expect_identical(mm_r1, c(1L, 6L))
})

test_that("complete landmark set grows leftward before righward", {
  mm_all <- landmarks_maxmin(bar, num = 6L, engine = "R")
  expect_identical(mm_all, c(1L, 6L, 3L, 2L, 4L, 5L))
})

# extended maxmin landmarks in R

test_that("'num' extensions extend numbers and preserve cover sets", {
  # starting with a number of landmarks
  mm_n <- landmarks_maxmin(bar, num = 2L, cover = TRUE)
  mm_n_nm <- landmarks_maxmin(bar, num = 2L, cover = TRUE,
                              extend_num = extension(mult = 1))
  mm_n_na <- landmarks_maxmin(bar, num = 2L, cover = TRUE,
                              extend_num = extension(add = 3))
  expect_true(all(mapply(identical,
                         x = mm_n$cover_set,
                         y = mm_n_nm$cover_set[seq_along(mm_n$cover_set)])))
  expect_true(all(mapply(identical,
                         x = mm_n_nm$cover_set,
                         y = mm_n_na$cover_set[seq_along(mm_n_nm$cover_set)])))
  # starting with a ball radius
  mm_r <- landmarks_maxmin(bar, radius = .4, cover = TRUE)
  mm_r_nm <- landmarks_maxmin(bar, radius = .4, cover = TRUE,
                              extend_num = extension(mult = .25))
  mm_r_na <- landmarks_maxmin(bar, radius = .4, cover = TRUE,
                              extend_num = extension(add = 3))
  expect_true(all(mapply(identical,
                         x = mm_r$cover_set,
                         y = mm_r_nm$cover_set[seq_along(mm_r$cover_set)])))
  expect_true(all(mapply(identical,
                         x = mm_r_nm$cover_set,
                         y = mm_r_na$cover_set[seq_along(mm_r_nm$cover_set)])))
})

test_that("'radius' extensions preserve numbers and extend cover sets", {
  # starting with a number of landmarks
  mm_n <- landmarks_maxmin(bar, num = 2L, cover = TRUE)
  mm_n_rm <- landmarks_maxmin(bar, num = 2L, cover = TRUE,
                              extend_radius = extension(mult = 1))
  mm_n_ra <- landmarks_maxmin(bar, num = 2L, cover = TRUE,
                              extend_radius = extension(add = 3))
  expect_true(all(mapply(function(x, y) all(x %in% y),
                         x = mm_n$cover_set,
                         y = mm_n_rm$cover_set)))
  expect_true(all(mapply(function(x, y) all(x %in% y),
                         x = mm_n_rm$cover_set,
                         y = mm_n_ra$cover_set)))
  # starting with a ball radius
  mm_r <- landmarks_maxmin(bar, radius = .4, cover = TRUE)
  mm_r_rm <- landmarks_maxmin(bar, radius = .4, cover = TRUE,
                              extend_radius = extension(mult = .25))
  mm_r_ra <- landmarks_maxmin(bar, radius = .4, cover = TRUE,
                              extend_radius = extension(add = 3))
  expect_true(all(mapply(function(x, y) all(x %in% y),
                         x = mm_r$cover_set,
                         y = mm_r_rm$cover_set)))
  expect_true(all(mapply(function(x, y) all(x %in% y),
                         x = mm_r_rm$cover_set,
                         y = mm_r_ra$cover_set)))
})

# firstlast

test_that("firstlast procedure prefers marginalizes multiple-points", {
  # both singleton points are centers because they are not the double-point
  expect_identical(firstlast(peg), c(1L, 2L))
  # only zero is a center because it is equi-nearest from the endpoints
  expect_identical(firstlast(bar), c(3L))
})

test_that("lastfirst procedure ...", {
  expect_identical(lastfirst(peg), c(3L, 4L))
  expect_identical(lastfirst(peg, ties_method = "max"), c(1L, 2L))
  # lastfirst set consists of points equidistant from no other pairs of points
  expect_identical(lastfirst(bar), c(1L, 4L, 6L))
  expect_identical(lastfirst(bar, ties_method = "max"), c(2L, 5L))
})

# -+- use 2-dimensional data sets -+-
test_that("firstlast procedure works on multiple data matrices", {
  expect_identical(firstlast(peg, bar + 1), c(1L))
  expect_identical(firstlast(peg, bar + 1.25), c(1L))
  expect_identical(firstlast(peg, bar + 1.5), c(2L))
  expect_identical(firstlast(peg, bar + 2), c(2L))
  expect_identical(firstlast(peg, bar + 4), c(3L, 4L))
})

# -+- use 2-dimensional data sets -+-
test_that("lastfirst procedure works on multiple data matrices", {
  expect_identical(lastfirst(peg, bar), seq(4))
  expect_identical(lastfirst(peg, bar + 1), c(2L, 3L, 4L))
  expect_identical(lastfirst(peg, bar + 2), c(1L, 3L, 4L))
  expect_identical(lastfirst(peg, bar + 3), seq(4))
  expect_identical(lastfirst(peg, bar + 4), c(1L, 2L))
})

# lastfirst landmarks in R

test_that("landmark set is obtained starting from duplicate point", {
  expect_error(landmarks_lastfirst(peg, seed_index = 3L), NA)
})

test_that("invalid numbers of sets prompt warnings", {
  expect_warning(landmarks_lastfirst(peg, num = 1L, cardinality = 2L),
                 "cardinality")
  expect_warning(landmarks_lastfirst(peg, num = 4L),
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
  lf_all <- suppressWarnings(landmarks_lastfirst(peg, num = 4L,
                                                 seed_index = 1L))
  expect_identical(lf_all, c(1L, 3L, 2L))
})

test_that("2-set cover of bar uses endpoints", {
  lf_n2 <- landmarks_lastfirst(bar, num = 2L)
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
  lf_all <- landmarks_lastfirst(bar, num = 6L)
  expect_identical(lf_all, c(1L, 6L, 3L, 4L, 2L, 5L))
})

# extended lastfirst landmarks in R

test_that("'num' extensions extend numbers and preserve cover sets", {
  # starting with a number of landmarks
  mm_n <- landmarks_lastfirst(bar, num = 2L, cover = TRUE)
  mm_n_nm <- landmarks_lastfirst(bar, num = 2L, cover = TRUE,
                                 extend_num = extension(mult = 1))
  mm_n_na <- landmarks_lastfirst(bar, num = 2L, cover = TRUE,
                                 extend_num = extension(add = 3))
  expect_true(all(mapply(identical,
                         x = mm_n$cover_set,
                         y = mm_n_nm$cover_set[seq_along(mm_n$cover_set)])))
  expect_true(all(mapply(identical,
                         x = mm_n_nm$cover_set,
                         y = mm_n_na$cover_set[seq_along(mm_n_nm$cover_set)])))
  # starting with a neighborhood cardinality
  mm_k <- landmarks_lastfirst(bar, cardinality = 2L, cover = TRUE)
  mm_k_nm <- landmarks_lastfirst(bar, cardinality = 2L, cover = TRUE,
                                 extend_num = extension(mult = .25))
  mm_k_na <- landmarks_lastfirst(bar, cardinality = 2L, cover = TRUE,
                                 extend_num = extension(add = 3))
  expect_true(all(mapply(identical,
                         x = mm_k$cover_set,
                         y = mm_k_nm$cover_set[seq_along(mm_k$cover_set)])))
  expect_true(all(mapply(identical,
                         x = mm_k_nm$cover_set,
                         y = mm_k_na$cover_set[seq_along(mm_k_nm$cover_set)])))
})

test_that("'cardinality' extensions preserve numbers and extend cover sets", {
  # starting with a number of landmarks
  mm_n <- landmarks_lastfirst(bar, num = 2L, cover = TRUE)
  mm_n_ka <- landmarks_lastfirst(bar, num = 2L, cover = TRUE,
                                 extend_cardinality = extension(add = 2))
  mm_n_km <- landmarks_lastfirst(bar, num = 2L, cover = TRUE,
                                 extend_cardinality = extension(mult = 1))
  expect_true(all(mapply(function(x, y) all(x %in% y),
                         x = mm_n$cover_set,
                         y = mm_n_ka$cover_set)))
  expect_true(all(mapply(function(x, y) all(x %in% y),
                         x = mm_n_ka$cover_set,
                         y = mm_n_km$cover_set)))
  # starting with a neighborhood cardinality
  mm_k <- landmarks_lastfirst(bar, cardinality = 2L, cover = TRUE)
  mm_k_km <- landmarks_lastfirst(bar, cardinality = 2L, cover = TRUE,
                                 extend_cardinality = extension(mult = 1))
  mm_k_ka <- landmarks_lastfirst(bar, cardinality = 2L, cover = TRUE,
                                 extend_cardinality = extension(add = 3))
  expect_true(all(mapply(function(x, y) all(x %in% y),
                         x = mm_k$cover_set,
                         y = mm_k_km$cover_set)))
  expect_true(all(mapply(function(x, y) all(x %in% y),
                         x = mm_k_km$cover_set,
                         y = mm_k_ka$cover_set)))
})
