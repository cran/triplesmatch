test_that("aberrant ranks make sense on continuous data", {
  set.seed(200)
  ymat <- matrix(rnorm(30), nrow = 10)
  treated1 <- rep(c(TRUE, FALSE), 5)
  ab <- aberrantscores(ymat, 0.5, cutoff_dir = "less", tau = 0, treated1 = NULL)
  expect_equal(which.min(ymat), which.max(ab))
  expect_equal(sum(ab != 0), sum(ymat < 0.5))
})

test_that("aberrant ranks make sense on discrete data and with null effect nonzero", {
  set.seed(21)
  ymat <- matrix(sample(1:10, 30, replace = TRUE), nrow = 10)
  treated1 <- rep(c(TRUE, FALSE), 5)
  ab <- aberrantscores(ymat, 7, cutoff_dir = "greater", tau = 0, treated1 = NULL)
  expect_equal(which.max(ymat), which.max(ab))
  expect_equal(sum(ab != 0), sum(ymat > 7))

  ab1 <- aberrantscores(ymat, 7, cutoff_dir = "greater", tau = 2, treated1 = treated1)
  expect_true(sum(ab1 != 0) < sum(ymat > 7))
  # Look at the rows where treated units are in columns 2 and 3
  # If it used to be <= 9, it should now not be aberrant after the adjustment for tau
  expect_true(all(ab1[!treated1, 2:3][ymat[!treated1, 2:3] <= 9] == 0))
  # If it used to be 10, it should still be aberrant
  expect_true(all(ab1[!treated1, 2:3][ymat[!treated1, 2:3] == 9] >= 0))

  # Look at the rows where controls are in columns 2 and 3
  # Controls who used to be 8 should still be aberrant
  expect_true(all(ab1[treated1, 2:3][ymat[!treated1, 2:3] == 8] >= 0))

})
