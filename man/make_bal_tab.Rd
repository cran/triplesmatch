% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_bal_tab.R
\name{make_bal_tab}
\alias{make_bal_tab}
\title{Make covariate balance table}
\usage{
make_bal_tab(X, z, m, cov_names)
}
\arguments{
\item{X}{Covariate matrix}

\item{z}{Vector of treatment indicators. Must be in same order as `y`}

\item{m}{`m` element of the list returned from `triples()` function containing information
about matched individuals}

\item{cov_names}{Row names to use instead of the column names of X when returning the table}
}
\value{
Table displaying means of the treated and control groups before and after matching,
    as well as standardized differences before and after matching
}
\description{
Make covariate balance table
}
\examples{
# Generate some data
set.seed(8)
n <- 200
nt <- floor(n * 0.5)
nc <- n - nt
x <- c(rnorm(nt, 0, 1), rnorm(nc, 0.6, 1))
z <- c(rep(1, nt), rep(0, nc))
names(z) <- 1:length(z)
names(x) <- 1:length(x)
# Create some strata
ps <- glm(z ~ x, family = binomial)$fitted.values
ps_st <- cut(ps, c(0, quantile(ps, 1/3 * 1:2), 1), labels = 1:3)
# Create a distance matrix
dist <- dist_mahal(data.frame(x = x), z, ps_st)
# Construct the triples match
triplesm <- triples(cost = dist, z = z, st = ps_st, solver = "rlemon")
make_bal_tab(X = cbind(x, ps), z = z, m = triplesm$m, cov_names = c("x", "prop score"))
}
