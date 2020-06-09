context("integer values with multiplicity")

# artificial integer-valued data set with several duplicates

set.seed(0)
m <- matrix(sample(30L)[rnbinom(n = 60L, size = 5, prob = 1/3)])
print(table(m))

# maxmin landmarks

test_that("landmarks are generated", {
  # maxmin landmarks in C++
  expect_silent(landmarks_maxmin_cpp(m, num_sets = nrow(m)))
  expect_silent(landmark_maxmin(m, n = nrow(m)))
  # maxmin landmarks in R
  expect_silent(landmarks_maxmin(m, num_sets = nrow(m)))
  expect_silent(landmarks_maxmin_R(m, num_sets = nrow(m)))
})

test_that("landmark sets agree", {
  expect_equal(landmarks_maxmin_cpp(m, num_sets = nrow(m)),
               landmarks_maxmin_R(m, num_sets = nrow(m)))
})

# lastfirst landmarks

test_that("landmarks are generated", {
  # lastfirst landmarks in C++
  expect_silent(landmarks_lastfirst_cpp(m, num_sets = nrow(m)))
  # lastfirst landmarks in R
  expect_silent(landmarks_lastfirst_R(m, num_sets = nrow(m)))
})

test_that("landmark sets agree", {
  expect_equal(landmarks_lastfirst_cpp(m, num_sets = nrow(m)),
               landmarks_lastfirst_R(m, num_sets = nrow(m)))
})
