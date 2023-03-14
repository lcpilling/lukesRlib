# lukesRlib 0.1.5

* `tidy_ci()` now tidies up terms included as categorical variables. "as.factor(x_var)1" is replaced with "x_var-1". Can be disabled with `tidy_factors=FALSE`
* addition of "Plotting" category of functions, starting with `doCoxSplinePlot()` that nicely plots the output of a spline term in a CoxPH model

# lukesRlib 0.1.4.1

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
