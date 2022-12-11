context("sample from circle without multiplicity")

# artificial data set sampled uniformly from the circle

set.seed(0)
s <- tdaunif::sample_circle(n = 24L)
plot(s, asp = 1, pch = 1)

# maxmin landmarks

test_that("landmarks are generated", {
  # maxmin landmarks in C++
  expect_silent(landmarks_maxmin(s, num = nrow(s), engine = "C++"))
  expect_silent(landmark_maxmin(s, num = nrow(s)))
  # maxmin landmarks in R
  expect_silent(landmarks_maxmin(s, num = nrow(s), engine = "original"))
  expect_silent(landmarks_maxmin(s, num = nrow(s), engine = "R"))
})

test_that("landmark sets agree", {
  expect_equal(landmarks_maxmin(s, num = nrow(s), engine = "original"),
               landmarks_maxmin(s, num = nrow(s), engine = "C++"),
               landmarks_maxmin(s, num = nrow(s), engine = "R"))
})

test_that("landmark sets and cover sets agree", {
  expect_equal(
    landmarks_maxmin(s, num = nrow(s), engine = "C++", cover = TRUE),
    landmarks_maxmin(s, num = nrow(s), engine = "R", cover = TRUE)
  )
})

# lastfirst landmarks

test_that("landmarks are generated", {
  # lastfirst landmarks in C++
  expect_silent(landmarks_lastfirst(s, num = nrow(s), engine = "C++"))
  # lastfirst landmarks in R
  expect_silent(landmarks_lastfirst(s, num = nrow(s), engine = "R"))
})

test_that("landmark sets and cover sets agree", {
  expect_equal(
    landmarks_lastfirst(s, num = nrow(s), engine = "C++", cover = TRUE),
    landmarks_lastfirst(s, num = nrow(s), engine = "R", cover = TRUE)
  )
})
