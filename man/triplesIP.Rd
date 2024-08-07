% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/triplesIP.R
\name{triplesIP}
\alias{triplesIP}
\title{Integer program for finding optimal triples match}
\usage{
triplesIP(z, cost, mt, mc, time_limit = Inf, threads = 1, verbose = 0)
}
\arguments{
\item{z}{Treatment indicator vector. 0 for control, 1 for treated}

\item{cost}{Matrix of costs. Rows correspond to treated units; columns to controls}

\item{mt}{The number of treated units to be used}

\item{mc}{The number of control units to be used}

\item{time_limit}{The amount of time in seconds before the solver should abort}

\item{threads}{The number of threads that should be allocated}

\item{verbose}{Whether the output of the 'gurobi' solver should be printed. 0 if not, 1 if so}
}
\value{
A named list with two elements: `match` and `opt_info`. `match` contains the triples match.
Similarly to the [triples()] function, this is in the form of a
 data.frame with number of rows equal to the number of triples and 8 columns specifying the
 match number, the names of the three units within the match, the costs of the two
 treated-control pairs within the match, the number of treated units, and the stratum.
 `opt_info` contains technical output from the optimization solver.
}
\description{
This finds the optimal triples match using a quadratic program. The 'gurobi'
package should be installed if using this function. This function should not be
used for large problems. Note that this solver may find a good solution even if not
optimal; setting `time_limit` is recommended. For most problems, [triples()] should
be used instead to find a good approximate solution very quickly.
}
\examples{
\dontshow{if (rlang::is_installed(c("gurobi"))) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Generate some data
set.seed(1)
n <- 40
x <- rnorm(n, 0, 1)
nt <- floor(n * 0.4)
nc <- n - nt
z <- c(rep(1, nt), rep(0, nc))
# Create some strata
ps <- glm(z ~ x, family = binomial)$fitted.values
ps_st <- cut(ps, c(0, quantile(ps, 1/3 * 1:2), 1), labels = 1:3)
# Create a distance matrix
dist <- dist_mahal(data.frame(x = x), z, ps_st)
# Construct the triples match using integer program for stratum 1
mIP <- triplesIP(z = z[ps_st == 1], cost = dist[[1]],
                 mt = 5, mc = 7, time_limit = 30, threads = 1, verbose = 0)
\dontshow{\}) # examplesIf}
}
\seealso{
triples for an approximate solution
}
