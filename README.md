<img align="right" src="https://github.com/lukepilling/lukesRlib/raw/master/lukesRlib.png" width="180" />

# lukesRlib
My library of R functions I sometimes find useful

[![](https://img.shields.io/badge/version-0.1.4-informational.svg)](https://github.com/lukepilling/lukesRlib)
[![](https://img.shields.io/github/last-commit/lukepilling/lukesRlib.svg)](https://github.com/lukepilling/lukesRlib/commits/master)
[![](https://img.shields.io/badge/lifecycle-experimental-9cf.svg)](https://www.tidyverse.org/lifecycle/#experimental)

<sub>Toolbox icon from https://vectorified.com/icon-tool-box</sub>

## List of functions
  - [Hypothesis testing](#hypothesis-testing)
    - [tidy_ci()](#tidy_ci)
    - [get_assoc()](#get_assoc)
  - [Data Transformation](#data-transformation)
    - [carrec()](#carrec)
    - [inv_norm()](#inv_norm)
    - [z_trans()](#z_trans)
  - [Working with test statistics](#working-with-test-statistics)
    - [get_se()](#get_se), [get_z()](#get_z), and [get_p()](#get_p)
    - [get_p_extreme()](#get_p_extreme), [get_p_neglog10()](#get_p_neglog10), [get_p_neglog10_n()](#get_p_neglog10_n)

## Installation
To install `lukesRlib` from GitHub use the `remotes` package:

`remotes::install_github("lukepilling/lukesRlib")`

To update the package just run the above command again.


## Hypothesis testing

### tidy_ci()
This function `tidy_ci()` runs [`broom::tidy()`](https://broom.tidymodels.org/) and returns the tidy estimates with CIs calculated as EST +/- 1.96*SE

Motivation: by default the [`broom`](https://broom.tidymodels.org/) package uses `confint()` to estimate CIs. For GLMs this calculates CIs via the profile likelihood method. When using large datasets this takes a long time and does not meaningfully alter the CIs compared to calculating using 1.96*SE

`tidy_ci()` does a few other nice things: hides the intercept by default, automatically detects logistic/CoxPH/CRR models and exponentiates the estimates, and if p==0 returns the 'extreme p' as a string. Other options include -log10 p-values.

See the [`tidy_ci()` Wiki page](https://github.com/lukepilling/lukesRlib/wiki/tidy_ci()) page for more details 

#### Examples

```R
fit_linear = glm(bmi ~ age + sex, data = d)
tidy_ci(fit_linear)
#> # A tibble: 2 x 8
#>   term  estimate std.error statistic   p.value conf.low conf.high p.extreme
#>   <chr>    <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl> <chr>    
#> 1 age     0.0196  0.000847     23.1  4.72e-118   0.0179    0.0212 NA       
#> 2 sex     0.703   0.0137       51.4  0           0.676     0.729  9.39e-576

library(survival)
fit_coxph = coxph(Surv(time, status) ~ age + sex + as.factor(smoking_status), data = d)
tidy_ci(fit_coxph, neglog10p=TRUE)
# CoxPH model :. estimate=exp()
#> # A tibble: 4 x 8
#>   term                       estimate std.error statistic  p.value conf.low conf.high neglog10p
#>   <chr>                         <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>     <dbl>
#> 1 age                           0.995  0.000837     -6.56 5.28e-11    0.993     0.996     10.3 
#> 2 sex                           1.04   0.0109        3.66 2.52e- 4    1.02      1.06       3.60
#> 3 as.factor(smoking_status)1    1.04   0.0120        3.26 1.13e- 3    1.02      1.06       2.95
#> 4 as.factor(smoking_status)2    1.03   0.0149        2.16 3.08e- 2    1.00      1.06       1.51

# ^^ automatically identified the input as from a coxph model and exponentiated estimate/CIs
```


### get_assoc()

To easily get tidy model output for a categorical or continuous exposure, including sample size (and cases ig logistic), outcome, and model info. Idea is to make quick loops easy. For all exposures, it gets the N. For categorical exposures, the N is split by group, and a row is included for the reference category. See the [`get_assoc()` Wiki page](https://github.com/lukepilling/lukesRlib/wiki/get_assoc()) page for details 

#### Examples

```R
# for one outcome, equivalent to `tidy_ci(glm(weight ~ height +age+sex, d=ukb))` - with added `n`
get_assoc(y="weight", x="height", z="+age+sex", d=ukb)
#  A tibble: 1 x 10
#>   outcome exposure estimate std.error statistic   p.value conf.low conf.high     n model
#>   <chr>   <chr>       <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl> <int> <chr>
#> 1 weight  height      0.762    0.0296      25.7 3.36e-137    0.704     0.820  4981 lm

# categorical exposure and startified analysis, with note
get_assoc(y="chd", x="smoking_status", z="+age", d=ukb |> filter(sex==1), logistic=TRUE, af=TRUE, note="Males only")
#  A tibble: 3 x 12
#>   outcome  exposure         estimate std.error statistic p.value conf.low conf.high     n n_cases model    note      
#>   <chr>    <chr>               <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl> <dbl>   <dbl> <chr>    <chr>     
#> 1 chd      smoking_status_0    NA       NA         NA    NA        NA         NA     1073     146 logistic Males only
#> 2 chd      smoking_status_1     1.24     0.126      1.72  0.0852    0.970      1.59   918     180 logistic Males only
#> 3 chd      smoking_status_2     1.48     0.181      2.16  0.0311    1.04       2.11   285      52 logistic Males only

# multiple exposures on single outcome, then combine output (i.e., a "PheWAS")
x_vars = c("bmi","ldl","sbp_0_avg")
res = do.call(rbind, lapply(x_vars, get_assoc, y="chd", z="+age+sex", d=ukb, logistic=TRUE))
#> # A tibble: 3 x 11
#>   outcome exposure  estimate std.error statistic  p.value conf.low conf.high     n n_cases model   
#>   <chr>   <chr>        <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl> <int>   <int> <chr>   
#> 1 chd     bmi          1.08    0.0113      6.97  3.10e-12    1.06       1.11  4692     324 logistic
#> 2 chd     ldl          0.980   0.0689     -0.300 7.64e- 1    0.856      1.12  4498     311 logistic
#> 3 chd     sbp_0_avg    1.01    0.00333     3.53  4.23e- 4    1.01       1.02  4561     316 logistic
```


## Data Transformation

### carrec()

From Steve Miller's package https://github.com/svmiller/stevemisc

`carrec()` (phonetically: "car-wreck") is a simple port of `car::recode()` to avoid clashes in the `{car}` package. This package offers STATA-like recoding features for R.

For example, assume a variable of interest is on a 1-10 scale. You want to code values 6 and above to be 1, and code values of 1-5 to be 0. Here’s how you would do that.

``` r
x <- seq(1, 10)
x
#>  [1]  1  2  3  4  5  6  7  8  9 10

carrec(x, "1:5=0;6:10=1")
#>  [1] 0 0 0 0 0 1 1 1 1 1
```

### inv_norm()

Inverse (quantile) normalise a quantitative trait (vector) i.e., transform to a normal distribution with mean=0 and sd=1

```r
x_in = inv_norm(x)

df = df |> mutate(x_in = inv_norm(x))
```

### z_trans()

Z-transform a quantitative trait (vector) i.e., convert to mean=0 and sd=1, maintaining original distribution

```r
x_z = z_trans(x)

df = df |> mutate(x_z = z_trans(x))
```

## Working with test statistics

### get_se()

Calculate Standard Error from Confidence Intervals.  

```r
lci = 0.1
uci = 0.3

# Default denominator is (1.96*2) equivalent to 95% confidence (p<0.05)
get_se(lci, uci)
#>  [1] 0.05102041

# custom denominator e.g., if CIs correspond to a p-value 5*10-8
get_se(lci, uci, denominator=lukesRlib::get_z(5e-8)*2)   
#>  [1] 0.01834421
```


### get_z()

Return a Z-statistic from a given p-value

```r
p = 1e-10
get_z(p)
#>  [1] 6.466951
```


### get_p()

Return a p-value from a z (or t) statistic

```r
z = 10
get_p(z)
#>  [1] 1.523971e-23
```


### get_p_extreme()

Normally R will round numbers < 1*10-324 to zero. This function returns the "extreme p-value" as a string. Provide a z (or t) statistic

```r
z = 50
get_p_extreme(z)
#>  [1] "2.16e-545"
```


### get_p_neglog10()

Returns the -log10 p-value. Provide a z (or t) statistic

```r
z = 50
get_p_neglog10(z)
#>  [1] 544.0977
```


### get_p_neglog10_n()

Returns the -log10 p-value. Provide a z (or t) statistic and n (sample size)

```r
z = 50
n = 100000
get_p_neglog10_n(z, n)
#>  [1] 537.9851
```
