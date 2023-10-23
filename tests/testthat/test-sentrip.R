test_that("sentrip and senfm give same results for m-scores", {

  set.seed(246)
  n <- 30
  x <- rnorm(n, 0, 3)
  nt <- floor(n * 0.5)
  nc <- n - nt
  z <- c(rep(1, nt), rep(0, nc))
  dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
  triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")

  y <- 1:40
  names(y) <- sample(1:40)

  ywide <- formattrip(m = triplesm, y = y, type = "wide")

  scores <- sensitivityfull::mscoref(ywide$ymat, ywide$treated1)
  # Undo the reversals in score for two treated subjects since sentrip will do this
  scores[!ywide$treated1, ] <- -scores[!ywide$treated1, ]

  res1 <- sentrip(scores = scores, treated1 = ywide$treated1, gamma = 1.2, alternative = "greater")
  res2 <- sensitivityfull::senfm(y = ywide$ymat, treated1 = ywide$treated1, gamma = 1.2)

  expect_equal(res1, res2)
})
