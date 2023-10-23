test_that("triples gives expected result", {
  set.seed(1)
  n <- 15
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.4)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")

  expect_equal(1 * sum(triplesm$m$nOfTreated == 1) + 2 * sum(triplesm$m$nOfTreated == 2), nt)
  expect_equal(2 * sum(triplesm$m$nOfTreated == 1) + 1 * sum(triplesm$m$nOfTreated == 2), nc)
  expect_equal(sum(!is.na(unique(c(triplesm$m$itreated, triplesm$m$jcontrol, triplesm$m$kthird)))),
    nt + nc)
  expect_equal(triplesm$obj, 2.8814972)
  expect_equal(triplesm$bound, 2.54455335)
})


test_that("triples gives expected result when all triples are 1:2", {
  set.seed(1)
  n <- 15
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.2)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")

  expect_equal(1 * sum(triplesm$m$nOfTreated == 1) + 2 * sum(triplesm$m$nOfTreated == 2), nt)
  expect_equal(2 * sum(triplesm$m$nOfTreated == 1) + 1 * sum(triplesm$m$nOfTreated == 2), 2 * nt)
  expect_equal(sum(!is.na(unique(c(triplesm$m$itreated, triplesm$m$jcontrol, triplesm$m$kthird)))),
               3 * nt)
  expect_equal(triplesm$obj, 1.1566469)
  expect_equal(triplesm$bound, triplesm$obj)
})



test_that("triples gives expected result when all triples are 2:1", {
  set.seed(10)
  n <- 15
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.8)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  expect_warning({triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")}, "6 treated will not be matched")

  expect_equal(1 * sum(triplesm$nOfTreated == 1) + 2 * sum(triplesm$m$nOfTreated == 2), min(nt, 2* nc))
  expect_equal(2 * sum(triplesm$m$nOfTreated == 1) + 1 * sum(triplesm$m$nOfTreated == 2), min(nc, 2 * nt))
  expect_equal(sum(!is.na(unique(c(triplesm$m$itreated, triplesm$m$jcontrol, triplesm$m$kthird)))), min(nt, 2*nc) + min(nc, 2*nt))
  expect_equal(triplesm$obj, 0.549446509)
  expect_equal(triplesm$obj, triplesm$bound)

})


test_that("triples gives expected result and bound when all triples are 2:1 but 1T and 1C must be dropped", {
  set.seed(10)
  n <- 20
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.65)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  expect_warning({triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")}, "1 treated individual will not be")
  t_to_use <- min(nt, 2* nc, 2 * floor((nc + nt)/3))

  expect_equal(1 * sum(triplesm$m$nOfTreated == 1) + 2 * sum(triplesm$m$nOfTreated == 2), t_to_use)
  expect_equal(2 * sum(triplesm$m$nOfTreated == 1) + 1 * sum(triplesm$m$nOfTreated == 2), t_to_use/2 )
  expect_equal(sum(!is.na(unique(c(triplesm$m$itreated, triplesm$m$jcontrol, triplesm$m$kthird)))),
               t_to_use * 3/2)
  expect_equal(triplesm$obj, 4.0939735)
  expect_equal(triplesm$bound, 2.62852155)

})



test_that("triples gives expected result when all triples are 2:1 but some controls must be dropped", {
  set.seed(10)
  n <- 20
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.62)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")
  t_to_use <- min(nt, 2* nc, 2 * floor((nc + nt)/3))

  expect_equal(1 * sum(triplesm$m$nOfTreated == 1) +
                 2 * sum(triplesm$m$nOfTreated == 2), t_to_use)
  expect_equal(2 * sum(triplesm$m$nOfTreated == 1) +
                 1 * sum(triplesm$m$nOfTreated == 2), t_to_use/2 )
  expect_equal(sum(!is.na(unique(c(triplesm$m$itreated, triplesm$m$jcontrol,
                                   triplesm$m$kthird)))), t_to_use * 3/2)
  expect_equal(triplesm$obj, 12.3374705)
  expect_true(triplesm$bound < triplesm$obj)
  expect_equal(triplesm$bound, 2.914628)
})

