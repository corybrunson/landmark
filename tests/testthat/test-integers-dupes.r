context("integer values with multiplicity")

# artificial integer-valued data set with several duplicates

set.seed(0)
multipoints <- matrix(sample(30L)[rnbinom(n = 60L, size = 5, prob = 1/3)])
print(table(multipoints))

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
