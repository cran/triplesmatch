% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infsentrip.R
\name{infsentrip}
\alias{infsentrip}
\title{Sensitivity analysis for triples matches informed by tests for unmeasured bias}
\usage{
infsentrip(gamma, sc, z, mset, alpha, alternative = "both")
}
\arguments{
\item{gamma}{The sensitivity parameter \eqn{\Gamma \geq 1}. Setting \eqn{\Gamma = 1}
performs a randomization test that assumes ignorable treatment assignment
given the matched triples}

\item{sc}{A matrix with N rows and at least two columns. The first column is the primary outcome
and the remaining columns are unaffected outcomes used to test for bias}

\item{z}{A vector of length N with 0s for control and 1s for treated units}

\item{mset}{A vector of length N indicating the matched triple}

\item{alpha}{A vector with length equal to the number of columns of `sc`. Each coordinate
contains the level of the test applied to the corresponding column of `sc`. If `alpha`
is a scalar, it is repeated for each column}

\item{alternative}{One of `greater`, `less` or `both`. `greater` implies the
alternative hypothesis is that the treatment has a positive effect on the scores of the primary outcome,
`less` implies the alternative hypothesis is that the treatment has a negative
effect on the scores of the primary outcome, and `both` conducts a two-sided hypothesis test.
The negative outcomes will always be two-sided tests (since one does not expect an effect in either direction)}
}
\value{
\describe{
\item{result }{Text indicating whether or not the test for bias rejects all biases of magnitude Gamma or less.  If yes, then the conclusion is that you must increase Gamma to continue.  If no, then the test on the primary outcome is conducted inside the confidence set defined by a test for bias.}
\item{optimization.problem }{Reiterates the result above, where the word yes means the optimization problem is infeasible, and the word no means it is feasible.  See the conclusion for a scientific interpretation of this aspect of the output.}
\item{conclusion }{Text indicating the result of the test for effect on the primary outcome.}
\item{deviates }{A vector of standardized deviates that might be compared with the standard Normal distribution.  There is one deviate for each column of `sc`.  If `sc` has column names, then the column names label the
  deviates.  The deviates are computed at the treatment assignment probabilities, theta, that solve the constrained optimization problem.}
\item{alphas }{A vector of significance levels used for the deviates, together with their total.  The total is relevant if the Bonferroni inequality is used to ensure joint level of all the tests.  The absolute deviates might be compared with qnorm(1-alphas/2) for a two-sided test.}
}
}
\description{
This function is very similar to [informedSen::informedsen()], with a few minor differences.
This version allows for matches to contain either two treated units and one control unit
or two controls and one treated unit
([informedSen::informedsen()] requires only one treated unit and a fixed number of controls in each match).
This version also allows the primary hypothesis test to be one-sided.
To use this function, the optimization software 'gurobi' and its R package must be installed.
}
\examples{
\dontshow{if (rlang::is_installed(c("gurobi", "Matrix"))) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Generate some data
set.seed(316)
n <- 30
x <- rnorm(n, 0, 1)
nt <- floor(n * 0.2)
nc <- n - nt
z <- c(rep(1, nt), rep(0, nc))
# Create a distance matrix (all units in one stratum here)
dist <- dist_mahal(data.frame(x = x), z, rep(1, n))[[1]]
# Conduct the triples match
triplesm <- triples_st(cost = dist, z = z, solver = "rlemon")
# Create primary and negative outcomes with some random unit names
y <- cbind(rnorm(40), runif(40))
rownames(y) <- sample(1:40)
# Reformat the triples match
ylong <- formattrip(m = triplesm, y = y, type = "long")
# Score the outcomes, in this case aberrant ranks for primary outcome and
#      ranks for unaffected outcome
y[, 1] <- aberrantscoreslong(y[, 1], cutoff = 0.5, cutoff_dir = "greater")
y[, 2] <- rank(y[, 2])
# Run the informed sensitivity analysis at gamma of 1.5
inf1 <- infsentrip(gamma = 1.5, sc = ylong$y, z = ylong$z, ylong$mset,
                  alpha = 0.05, alternative = "both")
\dontshow{\}) # examplesIf}
}
\seealso{
formattrip for easily creating inputs to this function.
}
