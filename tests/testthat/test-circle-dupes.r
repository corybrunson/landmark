context("sample from circle with multiplicity")

# artificial data set sampled with duplicates from the circle

set.seed(0)
x <- tdaunif::sample_circle(n = 12L)
y <- rbind(x, x[sample(nrow(x), 12L, replace = TRUE), , drop = FALSE])
plot(y, asp = 1, pch = 16, col = "#00000033")
l <- nrow(unique(y))

# maxmin landmarks

test_that("full landmark sets are generated", {
  # maxmin landmarks in C++
  expect_silent(landmarks_maxmin_cpp(y, num_sets = l))
  # maxmin landmarks in R
  expect_silent(landmarks_maxmin_R(y, num_sets = l))
})

test_that("warnings are generated", {
  expect_warning(landmarks_maxmin_cpp(y, num_sets = nrow(y)))
  expect_warning(landmarks_maxmin_R(y, num_sets = nrow(y)))
})

test_that("full landmark sets sets agree", {
  expect_equal(landmarks_maxmin_cpp(y, num_sets = l),
               landmarks_maxmin_R(y, num_sets = l))
})

# lastfirst landmarks

test_that("full landmarks are generated", {
  # lastfirst landmarks in C++
  expect_silent(landmarks_lastfirst_cpp(y, num_sets = l))
  # lastfirst landmarks in R
  expect_silent(landmarks_lastfirst_R(y, num_sets = l))
})

test_that("warnings are generated", {
  # lastfirst landmarks in C++
  expect_warning(landmarks_lastfirst_cpp(y, num_sets = nrow(y)))
  # lastfirst landmarks in R
  expect_warning(landmarks_lastfirst_R(y, num_sets = nrow(y)))
})

test_that("full landmark sets agree", {
  expect_equal(landmarks_lastfirst_cpp(y, num_sets = l),
               landmarks_lastfirst_R(y, num_sets = l))
})
