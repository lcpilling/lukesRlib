#' Get tidy model output with fast CIs 
#'
#' @description By default the (amazing) `broom` package uses the `confint()` function to calculate CIs. 
#' For GLMs this calculates confidence intervals via profile likelihood by default. 
#' When using large datasets this takes a long time and does not meaningfully alter the CIs 
#' compared to simply calculating using 1.96*SE.
#' This function `tidy_ci()` runs `broom::tidy()` and returns the tidy estimates with CIs 
#' calculated as EST +/- 1.96*SE.
#'
#' The function also does a few other nice/useful things to the output: hides the intercept by 
#' default, automatically detects logistic/CoxPH/CRR models and exponentiates the estimates, 
#' and if p=0 returns the extreme p as a string. Other optional outputs include -log10 p-values.
#'
#' Not tested for models other than `glm()` and `survival::coxph()` where it seems to work very well and produces consistent CIs. Also works well for `cmprsk::crr()`
#'
#' @return Returns a tibble - summary statistics from a model 
#'
#' @author Luke Pilling
#'
#' @name tidy_ci
#'
#' @param x object containing model output to be tidied e.g., from a `glm()` or `survival::coxph()`
#' @param ci calculate CIs using 1.96*SE method (default=TRUE)
#' @param intercept Exclude intercept for tidier output (default=FALSE)
#' @param extreme_ps If p=0 then return "extreme p-values" as strings (default=TRUE)
#' @param neglog10p Provides negative log10 p-values (if input is class `glm` or `coxph` or `crr` -- user can provide sample size `n=#` to override) (default=FALSE)
#' @param exp exponentiate estimate and CIs -- also see `check_model` (default=FALSE)
#' @param check_model set `exp=TRUE` if `glm(family=binomial)` or `survival::coxph()` or `cmprsk::crr()` was performed (default=TRUE)
#' @param n the N for `neglog10p` is extracted automatically for `glm` or `coxph` objects - override here if required (default=NA)
#' @param ... Other `tidy()` options
#'
#' @examples
#' fit_linear = glm(bmi ~ age + sex + as.factor(smoking_status), data = d)
#' tidy_ci(fit_linear)
#' 
#' fit_logistic = glm(current_smoker_vs_never ~ age + sex + bmi, data = d, family = binomial(link="logit"))
#' tidy_ci(fit_logistic)   # detect model and exponentiate automatically
#' tidy_ci(fit_logistic, check_model=FALSE)  # override auto checking to get untransformed estimates
#' 
#' fit_coxph = coxph(Surv(time, status) ~ age + sex + as.factor(smoking_status), data = d)
#' tidy_ci(fit_coxph)
#'
#' @export
#'

tidy_ci = function(x, 
                   ci = TRUE, 
                   exp = FALSE, 
                   intercept = FALSE, 
                   extreme_ps = TRUE,
                   neglog10p = FALSE, 
                   check_model = TRUE,
                   n = NA, 
                   conf.int = FALSE,     ## tidy() option
                   ...) {
	
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
			ret = ret |> dplyr::mutate(neglog10p=lukesRlib::get_p_neglog10_n(statistic, !!n))
			cat(paste0("N=", n, "\n"))
		}
	}
	
	## get extreme p-values?
	if (extreme_ps)  {
		if (any(ret$p.value==0))  {
			ret = ret |> mutate(p.extreme=if_else(p.value==0, lukesRlib::get_p_extreme(statistic), NA_character_))
		}
	}
	
	## exponentiate estimate and CIs?
	if (check_model & !exp)  {
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
