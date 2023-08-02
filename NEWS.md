# lukesRlib 0.2.5

* `get_assoc()` now expects the data as first argument, if not defined. Brings it more into the tidyverse. For instance, one can now:
  * `data |> filter(age>60) |> get_assoc(x="BMI", y="diabetes", z="sex")`
* Add theme `theme_minimal_modified` for ggplots
* Moves most dependencies to "Imports" - this is the preferred way to use existing packages https://r-pkgs.org/dependencies-mindset-background.html#sec-dependencies-namespace

# lukesRlib 0.2.4

* Add option "get_fit" to `get_assoc()` to also get model fit statistic when PheWASing. For each model type the `fit_stat` is:
  * lm = R-squared
  * logistic = McFadden's pseudo-R2
  * coxph = Harrell's c-statistic
* Fix error if user does not want to provide any covariates to `get_assoc()`
* `tidy_ci()` now returns the C-statistic for CoxPH models by default

# lukesRlib 0.2.3

* Decided to remove genetics-related functions to instead have them in a separate package, as the dependencies were getting out of hand
  * See https://github.com/lukepilling/gwasRtools
* Add `annotate_textp()` function for use in ggplot2 plots for adding easy text annotation geoms
* Moved some functions to be internal only to tidy up the namespace
* Fix `get_p_extreme` if NaN is passed rather than NA

# lukesRlib 0.2.2

* Decided to remove `doCoxSplinePlot()` because it is so much easier to do in {ggplot2} that a whole function is no longer necessary
  * Therefore the {pspline} dependency can be removed
* Set `extreme_ps` default to FALSE in `get_assoc()`

# lukesRlib 0.2.1

* Added p-value threshold option to `get_loci()`
* Fix problem getting `n` in `get_assoc()` when `model="coxph"`

# lukesRlib 0.2.0

* `get_assoc()` can now perform `coxph()` models from the {survival} package. The user provides strings depicting the survival object as the outcome e.g., "Surv(time, event)"
* Added two genetics-related functions: `lambda_gc()` (to calculate Lambda GC!) and `get_loci()` (to crudely determine independent loci in GWAS summary stats)
* General improvements/clarifications throughout code and documentation

# lukesRlib 0.1.9

* Added dependency {rlang} >= v1.0.0 so that the `sym()` function can be used - i.e., using provided variable name strings in `dplyr::mutate()` etc
* `get_assoc()` now allows options `inv_norm_x` and `inv_norm_y`

# lukesRlib 0.1.8

* Changed `get_assoc()` so that x and y can be vectors - i.e., multiple exposures/outcomes can be provided. Then uses {purrr} function `map2()` to iterate over combinations
* - This package now depends on {purrr} >= v1.0.0  &  {stringr} >= v1.5.0
* - Added functions `xv()` and `yv()` to improve ease of using {purrr} function `map2()` with `get_assoc()` for PheWAS of multiple exposures and outcomes in one analysis

# lukesRlib 0.1.7

* Calculate R^2 in `tidy_ci()` if `glm()` is linear regression
* Tidy up `tidy_ci()` parameter documentation

# lukesRlib 0.1.6

* fixes to `get_assoc()` (when no events in exposure group, or when two factor variables included)
* add `scale` options for exposure and outcome to `get_assoc()`
* in `tidy_ci()` tidy variable names with "scale()" 

# lukesRlib 0.1.5

* `tidy_ci()` now tidies up terms included as categorical variables. "as.factor(x_var)1" is replaced with "x_var-1". Can be disabled with `tidy_factors=FALSE`
* addition of "Plotting" category of functions, starting with `doCoxSplinePlot()` that nicely plots the output of a spline term in a CoxPH model
* Minor fix to `tidy_ci()` to avoid breaks due to missing p-values in model summary statistics

# lukesRlib 0.1.4

* Add new function `get_assoc()` to easily run PheWAS-like analysis using `tidy_ci()`

# lukesRlib 0.1.3

* Updates to `tidy_ci()` 
* Add `get_se()` from CIs function -- very simple but adds to collection of related functions
* General improvements

# lukesRlib 0.1.2

* Updated method for -log10 p-value calculation. One just uses Z. Other uses Z and N.

# lukesRlib 0.1.1

* Added two new p-value related functions -- very simple but nice to have available

# lukesRlib 0.1.0

* Added a `NEWS.md` file to track changes to the package.
