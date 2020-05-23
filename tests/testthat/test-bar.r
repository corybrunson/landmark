context("interval with different densities close to endpoints")

# maxmin landmarks

test_that("2-set cover uses endpoints", {
  lf_n2 <- landmarks_maxmin(bar, n = 2L)
  expect_identical(lf_n2, c(1L, 6L))
})

test_that("radius-1 cover uses endpoints only", {
  lf_e1 <- landmarks_maxmin(bar, eps = 1)
  expect_identical(lf_e1, c(1L, 6L))
})

test_that("complete landmark set grows leftward before righward", {
  lf_all <- landmarks_maxmin(bar, n = 6L)
  expect_identical(lf_all, c(1L, 6L, 3L, 2L, 4L, 5L))
})
