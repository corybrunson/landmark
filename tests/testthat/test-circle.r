context("sample from circle without multiplicity")

# artificial data set sampled uniformly from the circle

set.seed(0)
s <- tdaunif::sample_circle(n = 24L)
plot(s, asp = 1, pch = 1)

# maxmin landmarks

test_that("landmarks are generated", {
  # maxmin landmarks in C++
  expect_silent(landmarks_maxmin_cpp(s, num_sets = nrow(s)))
  # maxmin landmarks in R
  expect_silent(landmarks_maxmin(s, n = nrow(s)))
})

test_that("landmark sets agree", {
  expect_equal(landmarks_maxmin_cpp(s, num_sets = nrow(s)),
               landmarks_maxmin(s, n = nrow(s)))
})

# lastfirst landmarks

test_that("landmarks are generated", {
  # lastfirst landmarks in C++
  expect_silent(landmarks_lastfirst_cpp(s, num_sets = nrow(s)))
  # lastfirst landmarks in R
  expect_silent(landmarks_lastfirst_cory(s, number = nrow(s)))
})

test_that("landmark sets agree", {
  expect_equal(landmarks_lastfirst_cpp(s, num_sets = nrow(s)),
               landmarks_lastfirst_cory(s, number = nrow(s)))
})
