#' Get tidy model output with fast CIs 
#'
#' @description Function to run `broom::tidy()` and calculate CIs.
#'
#' By default the (amazing) \pkg{broom} package uses the `confint()` function to calculate CIs. 
#' For GLMs this calculates confidence intervals via profile likelihood by default. 
#' When using large datasets this takes a long time and does not meaningfully alter the CIs 
#' compared to simply calculating using 1.96*SE.
#' This function `tidy_ci()` runs `broom::tidy()` and returns the tidy estimates with CIs 
#' calculated as EST +/- 1.96*SE.
#'
#' The function also does a few other nice/useful things to the output: hides the intercept by 
#' default, calculates -log10 p-values, and automatically detects logistic/CoxPH/CRR models and 
#' exponentiates the estimates.
#'
#' Not tested for models other than `glm()` and `survival::coxph()` where it seems to work very well and produces consistent CIs. Also works well for `cmprsk::crr()`
#'
#' v0.20230109
#'
#' @return Returns a tibble - summary statistics from a model 
#'
#' @author Luke Pilling
#'
#' @name tidy_ci
#'
#' @param ci {default=TRUE} calculate CIs using 1.96*SE method
#' @param intercept {default=FALSE} Exclude intercept for tidier output
#' @param neglog10p {default=TRUE} Provides negative log10 p-values (if input is class `glm` or `coxph` or `crr` -- user can provide sample size `n=#` to override)
#' @param exp {default=FALSE} exponentiate estimate and CIs -- also see `check_family`
#' @param check_family {default=TRUE} set `exp=TRUE` if `glm(family=binomial)` or `survival::coxph()` or `cmprsk::crr()` was performed
#' @param n {default=NA} the N for `neglog10p` is extracted automatically for `glm` or `coxph` objects - override here if required
#' @param ... Other `tidy()` options
#'
#' @examples
#' fit_linear = glm(bmi ~ age + sex + as.factor(smoking_status), data = d)
#' tidy_ci(fit_linear)
#' 
#' fit_logistic = glm(current_smoker_vs_never ~ age + sex + bmi, data = d, family = binomial(link="logit"))
#' tidy_ci(fit_logistic)   # detect model and exponentiate automatically
#' tidy_ci(fit_logistic, check_family=FALSE)  # override auto checking to get untransformed estimates
#' 
#' fit_coxph = coxph(Surv(time_to_event, diagnosis_bin) ~ age + sex + bmi + as.factor(smoking_status), data = d)
#' tidy_ci(fit_coxph)
#'
#' @export
#'

tidy_ci = function(x = stop("Provide a model fit object"), 
		   ci = TRUE, 
		   exp = FALSE, 
		   intercept = FALSE, 
		   neglog10p = TRUE, 
		   check_family = TRUE,
		   n = NA, 
		   conf.int = FALSE,     ## tidy() option
		   ...) {
	
	#require(dplyr)
	#require(broom)
	
	## use `tidy()` CI method?  Only if not using the 1.96*SE method
	if (ci) conf.int = FALSE

	## get tidy output -- do not use `broom` CIs or Exponentiate options by default
	ret = broom::tidy(x, conf.int = conf.int, exponentiate = FALSE, ...)
	
	## get CIs based on 1.96*SE?
	if (ci)  ret = ret |> dplyr::mutate(conf.low=estimate-(1.96*std.error), conf.high=estimate+(1.96*std.error))
	
	## get -log10 p-value?
	if (neglog10p)  {
		if (is.na(n) & "glm" %in% class(x))  n = length(x$y)
		if (is.na(n) & "coxph" %in% class(x))  n = x$n
		if (is.na(n) & "crr" %in% class(x))  n = x$n
		if (is.na(n) & "tidycrr" %in% class(x))  n = x$cmprsk$n
		if (is.na(n)) cat("To calculate -log10 p-values provide the sample size `n`\n")
		if (!is.na(n)) {
			ret = ret |> dplyr::mutate(neglog10p=-1*(pt(abs(estimate/std.error),df=!!n,lower.tail=F,log.p=T) + log(2))/log(10))
			cat(paste0("N=", n, "\n"))
		}
	}
	
	## exponentiate estimate and CIs?
	if (check_family & !exp)  {
		model = ""
		if ("glm" %in% class(x)) {
			if (x$family$family == "binomial") {
				exp = TRUE
				model = "Binomial"
			}
		}
		if (any(c("coxph") %in% class(x)))  {
			exp = TRUE
			model = "CoxPH"
		}
		if (any(c("crr","tidycrr") %in% class(x)))  {
			exp = TRUE
			model = "CRR"
		}
		if (model != "")  cat(paste0(model, " model :: estimate=exp()\n"))
	}
	if (exp) ret = ret |> dplyr::mutate(estimate=exp(estimate))
	if (exp & (ci | conf.int)) ret = ret |> dplyr::mutate(conf.low=exp(conf.low), conf.high=exp(conf.high))
	
	## exclude intercept?
	if (!intercept) ret = ret |> dplyr::filter(term!="(Intercept)")
	
	## return object
	ret
	
}
