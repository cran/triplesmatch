test_that("long version of aberrant ranks is equivalent to wide version", {

  set.seed(21)
  y <- sample(1:10, 30, replace = TRUE)
  z <- c(rep(c(1, 0), 5), rep(c(0, 1), 10))
  ab_long <- aberrantscoreslong(y, 7, cutoff_dir = "greater", tau = 0, z = NULL)

  ymat <- matrix(y, nrow = 10)
  treated1 <- rep(c(TRUE, FALSE), 5)
  ab_wide <- aberrantscores(ymat, 7, cutoff_dir = "greater", tau = 0, treated1 = NULL)

  expect_equal(c(ab_wide), ab_long)

  ab1_long <- aberrantscoreslong(y, 7, cutoff_dir = "greater", tau = 2, z = z)
  ab1_wide <- aberrantscores(ymat, 7, cutoff_dir = "greater", tau = 2, treated1 = treated1)

  expect_equal(c(ab1_wide), ab1_long)
})
