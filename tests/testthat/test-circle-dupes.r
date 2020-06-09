context("sample from circle with multiplicity")

# artificial data set sampled with duplicates from the circle

set.seed(0)
x <- tdaunif::sample_circle(n = 12L)
m <- rbind(x, x[sample(nrow(x), 12L, replace = TRUE), , drop = FALSE])
plot(m, asp = 1, pch = 16, col = "#00000033")

# maxmin landmarks

test_that("landmarks are generated", {
  # maxmin landmarks in C++
  expect_silent(landmarks_maxmin_cpp(m, num_sets = nrow(m)))
  # maxmin landmarks in R
  expect_silent(landmarks_maxmin(m, n = nrow(m)))
})

test_that("landmark sets agree", {
  expect_equal(landmarks_maxmin_cpp(m, num_sets = nrow(m)),
               landmarks_maxmin(m, n = nrow(m)))
})

# lastfirst landmarks

test_that("landmarks are generated", {
  # lastfirst landmarks in C++
  expect_silent(landmarks_lastfirst_cpp(m, num_sets = nrow(m)))
  # lastfirst landmarks in R
  expect_silent(landmarks_lastfirst_cory(m, number = nrow(m)))
})

test_that("landmark sets agree", {
  expect_equal(landmarks_lastfirst_cpp(m, num_sets = nrow(m)),
               landmarks_lastfirst_cory(m, number = nrow(m)))
})
