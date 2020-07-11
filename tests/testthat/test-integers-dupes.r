context("integer values with multiplicity")

# artificial integer-valued data set with several duplicates

set.seed(0)
m <- matrix(sample(30L)[rnbinom(n = 60L, size = 5, prob = 1/3)])
print(table(m))
l <- length(unique(m))

# maxmin landmarks

test_that("full landmark sets are generated", {
  # maxmin landmarks in C++
  expect_silent(landmarks_maxmin(m, num_sets = l, engine = "C++"))
  # maxmin landmarks in R
  expect_silent(landmarks_maxmin(m, num_sets = l, engine = "R"))
})

test_that("warnings are generated", {
  expect_warning(landmarks_maxmin(m, num_sets = nrow(m), engine = "C++"))
  expect_warning(landmarks_maxmin(m, num_sets = nrow(m), engine = "R"))
})

test_that("full landmark sets and cover sets agree", {
  expect_equal(landmarks_maxmin(m, num_sets = l, engine = "C++", cover = TRUE),
               landmarks_maxmin(m, num_sets = l, engine = "R", cover = TRUE))
})

# lastfirst landmarks

test_that("full landmarks are generated", {
  # lastfirst landmarks in C++
  expect_silent(landmarks_lastfirst(m, num_sets = l, engine = "C++"))
  # lastfirst landmarks in R
  expect_silent(landmarks_lastfirst(m, num_sets = l, engine = "R"))
})

test_that("warnings are generated", {
  # lastfirst landmarks in C++
  expect_warning(landmarks_lastfirst(m, num_sets = nrow(m), engine = "C++"))
  # lastfirst landmarks in R
  expect_warning(landmarks_lastfirst(m, num_sets = nrow(m), engine = "R"))
})

test_that("full landmark sets and cover sets agree", {
  expect_equal(
    landmarks_lastfirst(m, num_sets = l, engine = "C++", cover = TRUE),
    landmarks_lastfirst(m, num_sets = l, engine = "R", cover = TRUE)
  )
})
