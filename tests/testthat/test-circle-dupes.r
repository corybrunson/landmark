context("sample from circle with multiplicity")

# artificial data set sampled with duplicates from the circle

set.seed(0)
X <- tdaunif::sample_circle(n = 12L)
multipoints <- rbind(X, X[sample(nrow(X), 12L, replace = TRUE), , drop = FALSE])
plot(multipoints, asp = 1, pch = 16, col = "#00000033")

# maxmin landmarks

test_that("landmarks are generated", {
  # maxmin landmarks in C++
  expect_silent(mm_lmk_cpp <- landmarks_maxmin_cpp(
    multipoints,
    num_sets = length(multipoints)
  ))
  # maxmin landmarks in R
  expect_silent(mm_lmk_r <- landmarks_maxmin(
    multipoints,
    n = length(multipoints)
  ))
})

test_that("landmark sets agree", {
  expect_equal(mm_lmk_cpp, mm_lmk_r)
})

# lastfirst landmarks

test_that("landmarks are generated", {
  # lastfirst landmarks in C++
  expect_silent(lf_lmk_cpp <- landmarks_lastfirst_cpp(
    multipoints,
    num_sets = length(multipoints)
  ))
  # lastfirst landmarks in R
  expect_silent(lf_lmk_r <- landmarks_lastfirst_cory(
    multipoints,
    number = length(multipoints)
  ))
})

test_that("landmark sets agree", {
  expect_equal(lf_lmk_cpp, lf_lmk_r)
})
