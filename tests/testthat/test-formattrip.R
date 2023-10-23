test_that("wide format established correctly", {

  set.seed(36)
  n <- 30
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.4)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")

  y <- 1:40
  names(y) <- sample(1:40)

  ywide <- formattrip(m = triplesm, y = y, type = "wide")

  # Everything is placed where we expect
  # Rownames are being used appropriately
  expect_true(all(y[as.character(triplesm$m$itreated[triplesm$m$nOfTreated == 1])] ==
                    ywide$ymat[triplesm$m$nOfTreated == 1, 1]))
  expect_true(all(y[as.character(triplesm$m$jcontrol[triplesm$m$nOfTreated == 1])] ==
                    ywide$ymat[triplesm$m$nOfTreated == 1, 2]))
  expect_true(all(y[as.character(triplesm$m$jcontrol[triplesm$m$nOfTreated == 2])] ==
                    ywide$ymat[triplesm$m$nOfTreated == 2, 1]))
  expect_true(all(y[as.character(triplesm$m$itreated[triplesm$m$nOfTreated == 2])] ==
                    ywide$ymat[triplesm$m$nOfTreated == 2, 2]))
  expect_true(all(y[as.character(triplesm$m$kthird)] ==
                    ywide$ymat[, 3]))
  expect_true(all(ywide$treated1[triplesm$m$nOfTreated == 1]))
  expect_false(all(ywide$treated1[triplesm$m$nOfTreated == 2]))

    })

test_that("error with no row names", {

  set.seed(36)
  n <- 30
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.4)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")

  y <- 1:40

  expect_error(formattrip(m = triplesm, y = y, type = "wide"), "names or rownames")
})

test_that("long format established correctly with one outcome", {

  set.seed(46)
  n <- 30
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.4)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")

  y <- 1:40
  names(y) <- sample(1:40)

  ylong <- formattrip(m = triplesm, y = y, type = "long")

  expect_true(all(ylong$y[ylong$mset == 1] %in%
        y[unlist(triplesm$m[triplesm$m$triple == 1, c("itreated", "jcontrol", "kthird")])]))

  two_treated <- which(triplesm$m$nOfTreated == 2)
  expect_true(all(ylong$y[ylong$mset == two_treated[1] & ylong$z == 1] %in%
                    y[unlist(triplesm$m[triplesm$m$triple == two_treated[1], c("itreated", "kthird")])]))
  expect_true(all(ylong$y[ylong$mset == two_treated[1] & ylong$z == 0] %in%
                    y[unlist(triplesm$m[triplesm$m$triple == two_treated[1], c("jcontrol")])]))
  })

test_that("error with no row names", {

  set.seed(46)
  n <- 30
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.4)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")

  y <- 1:40

  expect_error(formattrip(m = triplesm, y = y, type = "long"), "names or rownames")

  ymat <- cbind(y, rev(y))
  expect_error(formattrip(m = triplesm, y = ymat, type = "long"), "names or rownames")
})


test_that("long format established correctly with two outcomes", {
  set.seed(106)
  n <- 30
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.4)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")

  y <- cbind(1:40, runif(40))
  rownames(y) <- sample(1:40)

  ylong <- formattrip(m = triplesm, y = y, type = "long")

  expect_equal(ncol(ylong$y), 2)
  expect_true(all(ylong$y[ylong$mset == 1, ] %in%
                    y[unlist(triplesm$m[triplesm$m$triple == 1, c("itreated", "jcontrol", "kthird")]),]))

  two_treated <- which(triplesm$m$nOfTreated == 2)
  expect_true(all(ylong$y[ylong$mset == two_treated[1] & ylong$z == 1, ] %in%
                    y[unlist(triplesm$m[triplesm$m$triple == two_treated[1], c("itreated", "kthird")]), ]))
  expect_true(all(ylong$y[ylong$mset == two_treated[1] & ylong$z == 0, ] %in%
                    y[unlist(triplesm$m[triplesm$m$triple == two_treated[1], c("jcontrol")]), ]))
})

