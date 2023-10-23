test_that("triples gives expected result across strata", {
  set.seed(1)
  n <- 40
  x <- rnorm(n, 0, 1)
  nt <- floor(n * 0.4)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  ps <- glm(z ~ x, family = binomial)$fitted.values
  ps_st <- cut(ps, c(0, quantile(ps, 1/3 * 1:2), 1), labels = 1:3)

  dist <- dist_mahal(data.frame(x = x), z, ps_st)
  triplesm <- triples(cost = dist, z = z, st = ps_st, solver = "rlemon")
  c_drop <- sum(table(ps_st) %% 3)

  expect_equal(1 * sum(triplesm$m$nOfTreated == 1) + 2 * sum(triplesm$m$nOfTreated == 2), nt)
  expect_equal(2 * sum(triplesm$m$nOfTreated == 1) + 1 * sum(triplesm$m$nOfTreated == 2), nc - c_drop)
  expect_equal(sum(!is.na(unique(c(triplesm$m$itreated, triplesm$m$jcontrol,
                                   triplesm$m$kthird)))), nt + nc - c_drop)
  expect_equal(sum(triplesm$m$costStep1) + sum(triplesm$m$costStep2), 1.66256937)
})
