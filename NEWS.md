# lukesRlib 0.2.13 (21st Feb 2025)

### Bug fixes
* Fix tidy_ci() where concordance() from survival was not explicitly referenced properly
* Fix get_assoc() not properly returning the error when variables are missing


# lukesRlib 0.2.12 (16th Oct 2024)

### Bug fixes
* `get_assoc()` - Drop empty factor levels if using haven::as_factor()
* `get_assoc()` - Correctly identify covariate variable names when using haven::as_factor()

### Changes
* `get_assoc()` - New option "include_formula" allows user to request that the full regression formula (with covariates) is included in the output
* `get_assoc()` - Outputs exposure/outcome/covariate variables not in data.frame when `verbose=TRUE` (shows first 6 when FALSE)
* `get_assoc()` - Prints a bit more output to the screen by default


# lukesRlib 0.2.11 (30th Sept 2024)

### Changes
 * `tidy_ci()` now recognises and tidies factor naming if `haven::as_factor()` is used instead of the base R `as.factor()`
 * `get_assoc()` now uses `haven::as_factor()` instead of base R `as.factor()` for categorical exposures so that labels are shown nicely
 * Update dependencies to include {haven} and be more explicit which tidyverse packages are desired
 * Update example data categorical variables to include labels
 * Fix NEWS formatting


# lukesRlib 0.2.10 (6th June 2024)

### Changes
 * Change URL to reflect my GitHub username change from `lukepilling` to `lcpilling` to be more consistent between different logins, websites, and social media
   * https://lcpilling.github.io/lukesRlib
   * https://github.com/lcpilling/lukesRlib
 * `get_assoc()` 
   * Add options to Winsorize the exposure or outcome
   * Remove `beep`
   * Improved `verbose` reporting


# lukesRlib 0.2.9 (28th Feb 2024)

### Bug fixes
* `get_assoc()` - Fix `fit_stat_se` if missing 
* `get_assoc()` - Fix covariate checking if interaction included in formula. Resolving issue #1 


# lukesRlib 0.2.8 (19th Feb 2024)

* `get_assoc()` 
  * add `interacts_with` argument. Default is "". User can provide a variable name (string) to also get interaction terms
  * Fixed bug in `return_all_terms` when exposure was factor
  * Include SE for C-statistic if returning `fit_stat` for coxph() models

# lukesRlib 0.2.7 (29th Jan 2024)

* `get_assoc()` - add `return_all_terms` argument. Default is false. If TRUE, adds a new col `terms` and returns the estimates for all independent variables in the models.
* `get_assoc()` - remove `subset_d` argument
* `get_assoc()` - fixed if exposure or outcome needs to be protected by backticks (previously failed if variable was something like "bmi-1")
* `get_assoc()` - fixed problem where covariate estimates were also included in certain situations (e.g., if x="age" do not get covariate "percentage")
* Update examples for many functions

# lukesRlib 0.2.6 (13 Nov 2023)

* `get_assoc()` - fixed if no covariates provided
* `get_assoc()` - if x and y are the same, skip (useful if providing same list of variables twice for pairwise analysis)
* `get_assoc()` - added `verbose` option for bug checking
* `get_assoc()` - added `beep` option for beeping when finished
* `tidy_ci()` - add some functionality for summarizing fixed-effects from `lme4::lmer()`
* Add `example_data` object for examples. Contains ~5,000 participants from Framingham Heart Study

# lukesRlib 0.2.5 (25 Aug 2023)

* `get_assoc()` now expects the data as first argument, if not defined. Brings it more into the tidyverse. For instance, one can now:
  * `data |> filter(age>60) |> get_assoc(x="BMI", y="diabetes", z="sex")`
* Add theme `theme_minimal_modified` for ggplots
* Moves most dependencies to "Imports" - this is the preferred way to use existing packages https://r-pkgs.org/dependencies-mindset-background.html#sec-dependencies-namespace
* Fix -log10p calculation

# lukesRlib 0.2.4 (30 Jun 2023)

* Add option "get_fit" to `get_assoc()` to also get model fit statistic when PheWASing. For each model type the `fit_stat` is:
  * lm = R-squared
  * logistic = McFadden's pseudo-R2
  * coxph = Harrell's c-statistic
* Fix error if user does not want to provide any covariates to `get_assoc()`
* `tidy_ci()` now returns the C-statistic for CoxPH models by default

# lukesRlib 0.2.3 (7 Jun 2023)

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

# lukesRlib 0.2.0 (5 May 2023)

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

# lukesRlib 0.1.0 (26 Jan 2023)

* Added a `NEWS.md` file to track changes to the package.
