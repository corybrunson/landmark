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
  expect_silent(landmarks_maxmin(y, num = l, engine = "C++"))
  # maxmin landmarks in R
  expect_silent(landmarks_maxmin(y, num = l, engine = "R"))
})

test_that("warnings are generated", {
  expect_warning(landmarks_maxmin(y, num = nrow(y), engine = "C++"))
  expect_warning(landmarks_maxmin(y, num = nrow(y), engine = "R"))
})

test_that("full landmark sets and cover sets agree", {
  expect_equal(landmarks_maxmin(y, num = l, engine = "C++", cover = TRUE),
               landmarks_maxmin(y, num = l, engine = "R", cover = TRUE))
})

# lastfirst landmarks

test_that("full landmarks are generated", {
  # lastfirst landmarks in C++
  expect_silent(landmarks_lastfirst(y, num = l, engine = "C++"))
  # lastfirst landmarks in R
  expect_silent(landmarks_lastfirst(y, num = l, engine = "R"))
})

test_that("warnings are generated", {
  # lastfirst landmarks in C++
  expect_warning(landmarks_lastfirst(y, num = nrow(y), engine = "C++"))
  # lastfirst landmarks in R
  expect_warning(landmarks_lastfirst(y, num = nrow(y), engine = "R"))
})

test_that("full landmark sets and cover sets agree", {
  expect_equal(
    landmarks_lastfirst(y, num = l, engine = "C++", cover = TRUE),
    landmarks_lastfirst(y, num = l, engine = "R", cover = TRUE)
  )
})
