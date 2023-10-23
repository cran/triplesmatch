test_that("triplesIP gives same results as triples when all triples 2:1", {
  skip_if_not_installed("gurobi")
  set.seed(10)
  n <- 15
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.8)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  expect_warning({triplesm <- triples_st(cost = dist, z = z)}, "6 treated will not be matched")
  mIP <- triplesIP(z = z, cost = dist, mt = 6, mc = 3, time_limit = 60, threads = 1, verbose = 0)
  # Total cost matches the triples algo
  expect_equal(triplesm$obj, sum(mIP$match$costStep1) + sum(mIP$match$costStep2))
  # Total cost matches the stated objective value
  expect_equal(mIP$opt_info$objval, sum(mIP$match$costStep1) + sum(mIP$match$costStep2))
  # Stated objective value matches objective bound
  expect_equal(mIP$opt_info$objval, mIP$opt_info$objbound)
})

test_that("triplesIP gives same results as triples when all triples 1:2", {
  skip_if_not_installed("gurobi")
  set.seed(1)
  n <- 15
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.2)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z)

  mIP <- triplesIP(z = z, cost = dist, mt = 3, mc = 6, time_limit = 60, threads = 1, verbose = 0)
  # Total cost matches the triples algo
  expect_equal(triplesm$obj, sum(mIP$match$costStep1) + sum(mIP$match$costStep2))
  # Total cost matches the stated objective value
  expect_equal(mIP$opt_info$objval, sum(mIP$match$costStep1) + sum(mIP$match$costStep2))
  # Stated objective value matches objective bound
  expect_equal(mIP$opt_info$objval, mIP$opt_info$objbound)
})



test_that("triplesIP does better when it finishes", {
  skip_if_not_installed("gurobi")
  set.seed(6)
  n <- 30
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.4)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z)

  mt <- sum(triplesm$m$nOfTreated)
  mc <- nrow(triplesm$m) * 3 - mt
  mIP <- triplesIP(z = z, cost = dist, mt = mt, mc = mc, time_limit = 60, threads = 1, verbose = 0)
  # Objective beats the triples algo
  expect_true(mIP$opt_info$objval < triplesm$obj)
  # Stated objective value matches objective bound
  expect_equal(mIP$opt_info$objval, mIP$opt_info$objbound)
  # Objective exceeds lower bound of triples algo
  expect_true(mIP$opt_info$objval > triplesm$bound)
})
