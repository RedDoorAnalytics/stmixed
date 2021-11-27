# stmixed

## Multilevel mixed effects parametric survival analysis in Stata

The `stmixed` package began as the accompanying Stata implementation for the methodology in:

> Crowther MJ, Look MP, Riley RD. Multilevel mixed effects parametric survival models using adaptive Gauss-Hermite quadrature with application to recurrent events and individual participant data meta-analysis. *Statistics in Medicine* 2014;33(22):3844-3858.

`stmixed` has now become a wrapper function for `merlin`, which gives it a lot more capabilities. The new developments in `stmixed` are described in a recent *Stata Journal* paper:

> Crowther MJ. Multilevel mixed effects parametric survival analysis: Estimation, simulation and application. *Stata Journal* 2019;19(4):931-949. (Pre-print: https://arxiv.org/abs/1709.06633)

## Installation

The latest stable version of `stmixed` can be installed with:

```{stata}
ssc install stmixed
```

To install directly from this GitHub repository, use:

```{stata}
net install stmixed, from("https://raw.githubusercontent.com/RedDoorAnalytics/stmixed/main/")
```
