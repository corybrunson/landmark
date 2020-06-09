context("sample from circle with multiplicity")

# artificial data set sampled with duplicates from the circle

set.seed(0)
x <- tdaunif::sample_circle(n = 12L)
y <- rbind(x, x[sample(nrow(x), 12L, replace = TRUE), , drop = FALSE])
plot(y, asp = 1, pch = 16, col = "#00000033")

# maxmin landmarks

test_that("landmarks are generated", {
  # maxmin landmarks in C++
  expect_silent(landmarks_maxmin_cpp(y, num_sets = nrow(y)))
  # maxmin landmarks in R
  expect_silent(landmarks_maxmin(y, n = nrow(y)))
})

test_that("landmark sets agree", {
  expect_equal(landmarks_maxmin_cpp(y, num_sets = nrow(y)),
               landmarks_maxmin(y, n = nrow(y)))
})

# lastfirst landmarks

test_that("landmarks are generated", {
  # lastfirst landmarks in C++
  expect_silent(landmarks_lastfirst_cpp(y, num_sets = nrow(y)))
  # lastfirst landmarks in R
  expect_silent(landmarks_lastfirst_cory(y, number = nrow(y)))
})

test_that("landmark sets agree", {
  expect_equal(landmarks_lastfirst_cpp(y, num_sets = nrow(y)),
               landmarks_lastfirst_cory(y, number = nrow(y)))
})
