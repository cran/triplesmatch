test_that("infsentrip gives same results as informedsen if all triples 1:2 and two-sided", {
  skip_if_not_installed("gurobi")
  set.seed(316)
  n <- 30
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.2)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")

  y <- cbind(rnorm(40), runif(40))
  rownames(y) <- sample(1:40)

  ylong <- formattrip(m = triplesm, y = y, type = "long")

  inf1 <- infsentrip(gamma = 1.5, sc = ylong$y, z = ylong$z, ylong$mset, alpha = 0.05, alternative = "both")

  skip_if_not_installed("gurobi")
  skip_if_not_installed("Matrix")
  inf2 <- informedSen::informedsen(gamma = 1.5, sc = ylong$y, z = ylong$z, ylong$mset, alpha = 0.05)
  expect_equal(inf1, inf2, tolerance = 2e-3)
  })



test_that("infsentrip one sided is similar to informedsen", {
  skip_if_not_installed("gurobi")
  set.seed(336)
  n <- 30
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.2)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")

  y <- cbind(rnorm(40), runif(40))
  rownames(y) <- sample(1:40)

  ylong <- formattrip(m = triplesm, y = y, type = "long")
  ylong$y[ylong$z == 1, 1] <- ylong$y[ylong$z ==1, 1] + 0.25

  # One sided test rejects
  inf1 <- infsentrip(gamma = 1, sc = ylong$y, z = ylong$z, ylong$mset, alpha = 0.05, alternative = "greater")
  expect_match(inf1$conclusion, "rejects")

  # Doesnt reject if the wrong direction
  inf3 <- infsentrip(gamma = 1, sc = ylong$y, z = ylong$z, ylong$mset, alpha = 0.05, alternative = "less")
  expect_match(inf3$conclusion, "fails to reject")

})

