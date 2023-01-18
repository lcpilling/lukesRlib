# lukesRlib
Luke's library of R functions I sometimes find useful

[![](https://img.shields.io/badge/version-0.1.0-informational.svg)](https://github.com/lukepilling/lukesRlib)
[![](https://img.shields.io/github/last-commit/lukepilling/lukesRlib.svg)](https://github.com/lukepilling/lukesRlib/commits/master)
[![](https://img.shields.io/badge/lifecycle-experimental-9cf.svg)](https://www.tidyverse.org/lifecycle/#experimental)

## Table of Contents
  - [installation](#installation)
  - [tidy_ci()](#tidy_ci)
  - [carrec()](#carrec)
  - [inv_norm()](#inv_norm)
  - [z_trans()](#z_trans)
  - [get_extreme_p()](#get_extreme_p)
  - [get_neglog10_p()](#get_neglog10_p)

## installation
To install `lukesRlib` from GitHub use the `remotes` package:

`remotes::install_github("lukepilling/lukesRlib")`

To update the package just run the above command again.


## tidy_ci()
Function to run `broom::tidy()` and calculate CIs

By default the (amazing) `broom` package uses the `confint()` function to calculate CIs. For GLMs this calculates confidence intervals via profile likelihood by default. When using large datasets this takes a long time and does not meaningfully alter the CIs compared to simply calculating using 1.96*SE

This function `tidy_ci()` runs `broom::tidy()` and returns the tidy estimates with CIs calculated as EST +/- 1.96*SE

Also does a few other nice/useful things to the output: hides the intercept by default, calculates -log10 p-values, and automatically detects logistic/CoxPH/CRR models and exponentiates the estimates

### Options:
 - `ci` {default=TRUE} calculate CIs using 1.96*SE method
 - `intercept` {default=FALSE} Exclude intercept for tidier output
 - `neglog10p` {default=TRUE} Provides negative log10 p-values (if input is class `glm` or `coxph` or `crr` -- user can provide sample size `n=#` to override)
 - `exp` {default=FALSE} exponentiate estimate and CIs -- also see `check_family`
 - `check_family` {default=TRUE} set `exp=TRUE` if `glm(family=binomial)` or `survival::coxph()` or `cmprsk::crr()` was performed
 - `n` {default=NA} the N for `neglog10p` is extracted automatically for `glm` or `coxph` objects - override here if required
 - `...` Other `tidy()` options 

Not tested for models other than `glm()` and `survival::coxph()` where it seems to work very well and produces consistent CIs. Also works well for `cmprsk::crr()` (and therefore `tidycmprsk::crr()`)

### Examples

```R
fit_linear = glm(bmi ~ age + sex, data = d)
tidy_ci(fit_linear)
#> # A tibble: 4 x 8
#>   term                       estimate std.error statistic   p.value conf.low conf.high neglog10p
#>   <chr>                         <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>     <dbl>
#> 1 age                          0.0196  0.000847     23.1  4.72e-118   0.0179    0.0212     117. 
#> 2 sex                          0.703   0.0137       51.4  0           0.676     0.729      574. 

library(survival)
fit_coxph = coxph(Surv(time, status) ~ age + sex + as.factor(smoking_status), data = d)
tidy_ci(fit_coxph)
#> # A tibble: 5 x 8
#>   term                       estimate std.error statistic  p.value conf.low conf.high neglog10p
#>   <chr>                         <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>     <dbl>
#> 1 age                           0.995  0.000837     -6.56 5.28e-11    0.993     0.996     10.3 
#> 2 sex                           1.04   0.0109        3.66 2.52e- 4    1.02      1.06       3.60
#> 3 as.factor(smoking_status)1    1.04   0.0120        3.26 1.13e- 3    1.02      1.06       2.95
#> 4 as.factor(smoking_status)2    1.03   0.0149        2.16 3.08e- 2    1.00      1.06       1.51

# ^^ automatically identified the input as from a coxph model and exponentiated estimate/CIs
```

## carrec()

Stolen straight from Steve Miller's package https://github.com/svmiller/stevemisc

### `carrec()`: A Port of `car::recode()`

`carrec()` (phonetically: “car-wreck”) is a simple port of
`car::recode()` that I put in this package because of various function
clashes in the `{car}` package. For those who cut their teeth on Stata,
this package offers Stata-like recoding features that are tough to find
in the R programming language.

For example, assume the following vector that is some variable of
interest on a 1-10 scale. You want to code the variables that are 6 and
above to be 1 and code the variables of 1-5 to be 0. Here’s how you
would do that.

``` r
x <- seq(1, 10)
x
#>  [1]  1  2  3  4  5  6  7  8  9 10

carrec(x, "1:5=0;6:10=1")
#>  [1] 0 0 0 0 0 1 1 1 1 1
```

## inv_norm()

Inverse (quantile) normalise a quantitative trait (vector) i.e., transform to a normal distribution with mean=0 and sd=1

```r
x_in = inv_norm(x)

df = df |> mutate(x_in = inv_norm(x))
```

## z_trans()

Z-transform a quantitative trait (vector) i.e., convert to mean=0 and sd=1, maintaining original distribution

```r
x_z = z_trans(x)

df = df |> mutate(x_z = z_trans(x))
```

## get_extreme_p()

Return a p-value even if <1*10-324 (returns a string) -- provide a z (or t) statistic

```r
z = 50
get_extreme_p(z)
#>  [1] "p value is 2.16 times 10^(-545)"
```

## get_neglog10_p()

Returns the -log10 p-value. Provide a z (or t) statistic and n (sample size)

```r
z = 50
n = 100000
get_neglog10_p(z, n)
#>  [1] 537.9851
```
