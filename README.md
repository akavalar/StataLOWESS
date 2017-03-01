# StataLOWESS

## Cleveland's implementation of LOWESS for Stata

Stata implementation of LOWESS that follows the original Fortran code by Cleveland (making sure that boundary subsets are not decreasing in size). Original implementation in Stata uses subsets of decreasing size.

For details, see this [Cross Validated post](http://stats.stackexchange.com/questions/262720/lowess-implementation-in-stata-vs-r-and-python).

`ksm_original7.ado` is the `ksm.ado` file found on [Stata Corp's website](http://www.stata.com/support/updates/stata7/ado/).

## Illustration of the problem
![Illustration](/Example.png?raw=true "Series diverge in the tails")
