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
fit_linear = glm(bmi ~ age + sex + as.factor(smoking_status), data = d)
tidy_ci(fit_linear)
## A tibble: 4 x 8
#  term                       estimate std.error statistic   p.value conf.low conf.high neglog10p
#  <chr>                         <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>     <dbl>
#1 age                          0.0196  0.000847     23.1  4.72e-118   0.0179    0.0212     117. 
#2 sex                          0.703   0.0137       51.4  0           0.676     0.729      574. 
#3 as.factor(smoking_status)1   0.630   0.0149       42.3  0           0.601     0.659      390. 
#4 as.factor(smoking_status)2  -0.203   0.0228       -8.91 5.28e- 19  -0.248    -0.159       18.3

fit_logistic = glm(current_smoker_vs_never ~ age + sex + bmi, data = d, family = binomial(link="logit"))
tidy_ci(fit_logistic)   ## automatically identifies a logistic model and exponentiates estimate/CIs
## A tibble: 3 x 8
#  term  estimate std.error statistic   p.value conf.low conf.high neglog10p
#  <chr>    <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>     <dbl>
#1 age      0.983  0.000584    -29.2  7.19e-188    0.982     0.984     187. 
#2 sex      1.70   0.00962      55.4  0            1.67      1.74      665. 
#3 bmi      0.991  0.00103      -8.44 3.08e- 17    0.989     0.993      16.5

library(survival)
fit_coxph = coxph(Surv(time_to_event, diagnosis_bin) ~ age + sex + bmi + as.factor(smoking_status), data = d)
tidy_ci(fit_coxph)   ## automatically identifies a coxph model and exponentiates estimate/CIs
## A tibble: 5 x 8
#  term                       estimate std.error statistic  p.value conf.low conf.high neglog10p
#  <chr>                         <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>     <dbl>
#1 age                           0.995  0.000837     -6.56 5.28e-11    0.993     0.996     10.3 
#2 sex                           1.04   0.0109        3.66 2.52e- 4    1.02      1.06       3.60
#3 bmi                           0.994  0.00100      -6.50 7.92e-11    0.992     0.995     10.1 
#4 as.factor(smoking_status)1    1.04   0.0120        3.26 1.13e- 3    1.02      1.06       2.95
#5 as.factor(smoking_status)2    1.03   0.0149        2.16 3.08e- 2    1.00      1.06       1.51
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

## get_extreme_p()

Return a p-value even if <1*10-324 (returns a string) -- provide a z or statistic

```r
get_extreme_p(z)
```

## get_neglog10_p()

Returns the -log10 p-value. Provide a beta, se, and n (sample size)

```r
get_neglog10_p(beta, se, n)
```
