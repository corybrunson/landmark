context("sample from circle without multiplicity")

# artificial data set sampled uniformly from the circle

set.seed(0)
s <- tdaunif::sample_circle(n = 24L)
plot(s, asp = 1, pch = 1)

# maxmin landmarks

test_that("landmarks are generated", {
  # maxmin landmarks in C++
  expect_silent(mm_lmk_cpp <- landmarks_maxmin_cpp(s, num_sets = length(s)))
  # maxmin landmarks in R
  expect_silent(mm_lmk_r <- landmarks_maxmin(s, n = length(s)))
})

test_that("landmark sets agree", {
  expect_equal(mm_lmk_cpp, mm_lmk_r)
})

# lastfirst landmarks

test_that("landmarks are generated", {
  # lastfirst landmarks in C++
  expect_silent(lf_lmk_cpp <- landmarks_lastfirst_cpp(s, num_sets = length(s)))
  # lastfirst landmarks in R
  expect_silent(lf_lmk_r <- landmarks_lastfirst_cory(s, number = length(s)))
})

test_that("landmark sets agree", {
  expect_equal(lf_lmk_cpp, lf_lmk_r)
})
