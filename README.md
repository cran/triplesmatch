
# triplesmatch

<!-- badges: start -->
<!-- badges: end -->

The goal of triplesmatch is to attain excellent covariate balance by matching two treated units and one control unit or vice versa within strata. Using such triples, as
opposed to also allowing pairs of treated and control units, 
allows easier interpretation of the two possible 
weights of observations and better insensitivity to unmeasured bias in the test statistic. Using triples instead of matching in a fixed 1:2 or 2:1 ratio allows for the match to be feasible in more situations.

The rrelaxiv package, which provides an alternative solver for the underlying network flow problems, carries an
academic license and is not available on CRAN, but
may be downloaded from Github at 
<https://github.com/josherrickson/rrelaxiv/>.
The 'Gurobi' commercial optimization software is required to use the two functions ``infsentrip()` and `triplesIP()`. These functions are not essential to the function of this package.
 The 'gurobi' R package can be installed following the instructions at <https://www.gurobi.com/documentation/9.1/refman/ins_the_r_package.html>.

## Installation

You can install the release version of optrefine from [CRAN](https://cran.r-project.org/) with:

``` r
install.packages("triplesmatch")
```
