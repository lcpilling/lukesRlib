#' Get tidy model output with fast CIs 
#'
#' @description By default the (amazing) `broom` package uses the `confint()` function to calculate CIs. 
#' For GLMs this calculates confidence intervals via profile likelihood by default. 
#' When using large datasets this takes a long time and does not meaningfully alter the CIs 
#' compared to simply calculating using 1.96*SE.
#' This function `tidy_ci()` runs `broom::tidy()` and returns the tidy estimates with CIs 
#' calculated as EST +/- 1.96*SE. (Well, actually 1.959964 :: from `lukesRlib::get_z(0.05)`)
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
#' @param x object containing model output to be tidied e.g., from a `glm()` or `survival::coxph()`.
#' @param ci Logical. Default is TRUE. Calculate CIs using 1.96*SE method - where 1.96 can be modified using `ci_denominator`.
#' @param ci_denominator Numeric. Default is 1.96. The standard error of the sample mean (default actually 1.959964 :: from `get_z(0.05)`).
#' @param intercept Logical. Default is FALSE. Exclude intercept for tidier output.
#' @param tidy_factors Logical. Default is TRUE. Tidy `as.factor(x_var)#` terms to `x_var-#`.
#' @param extreme_ps Logical. Default is TRUE. If p=0 then return "extreme p-values" as strings.
#' @param neglog10p Logical. Default is FALSE. Provides negative log10 p-values (if input is class `glm` or `coxph` or `crr` -- user can provide sample size `n=#` to override).
#' @param exp Logical. Default is FALSE. Exponentiate estimate and CIs -- also see `check_model`.
#' @param check_model Logical. Default is TRUE. Set `exp=TRUE` if `glm(family=binomial)` or `survival::coxph()` or `cmprsk::crr()` was performed.
#' @param n Numeric. Default is NA. The N for `neglog10p` is extracted automatically for `glm` or `coxph` objects - override here if required.
#' @param get_r2 Logical. Default is TRUE. get R^2 value if model is linear.
#' @param quiet Logical. Default is FALSE. Suppress text output
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
                   ci_denominator = 1.959964,
                   exp = FALSE, 
                   intercept = FALSE, 
                   tidy_factors = TRUE,
                   extreme_ps = TRUE,
                   neglog10p = FALSE, 
                   check_model = TRUE,
                   n = NA, 
                   get_r2 = TRUE,
                   conf.int = FALSE,     ## tidy() option
                   quiet = FALSE,
                   ...) {
	
	# use `tidy()` CI method?  Only if not using the 1.96*SE method
	if (ci) conf.int = FALSE
	
	# get tidy output -- do not use `broom` CIs or Exponentiate options by default
	ret = broom::tidy(x, conf.int = conf.int, exponentiate = FALSE, ...)
	
	# exclude intercept?
	if (!intercept) ret = ret |> dplyr::filter(term!="(Intercept)")
	
	# tidy variable name if scaled
	ret = ret |> dplyr::mutate(
		term=if_else(grepl("scale(", term, fixed=TRUE), str_replace(term, fixed(")"), "-scaled"), term),
		term=str_replace(term, fixed("scale("), "")
	)
	
	# tidy factor names?
	if (tidy_factors)  {
		ret = ret |> dplyr::mutate(
			term=if_else(grepl("as.factor(", term, fixed=TRUE), str_replace(term, fixed(")"), "-"), term),
			term=str_replace(term, fixed("as.factor("), "")
		)
	}
	
	# get CIs based on 1.96*SE?
	if (ci)  ret = ret |> dplyr::mutate(conf.low=estimate-(!!ci_denominator*std.error), conf.high=estimate+(!!ci_denominator*std.error))
	
	# Print N?
	text_out = ""
	if (is.na(n) & "glm" %in% class(x))     n = length(x$y)
	if (is.na(n) & "coxph" %in% class(x))   n = x$n
	if (is.na(n) & "crr" %in% class(x))     n = x$n
	if (is.na(n) & "tidycrr" %in% class(x)) n = x$cmprsk$n
	if (!is.na(n)) text_out = paste0("N=", n)
	
	# get -log10 p-value?
	if (neglog10p)  {
		if (is.na(n)) cat("To calculate -log10 p-values provide the sample size `n`\n")
		if (!is.na(n)) ret = ret |> dplyr::mutate(neglog10p=lukesRlib::get_p_neglog10_n(statistic, !!n))
	}
	
	# get extreme p-values?
	if (extreme_ps)  if (any(ret$p.value==0, na.rm=TRUE))  ret = ret |> dplyr::mutate(p.extreme=dplyr::if_else(p.value==0, lukesRlib::get_p_extreme(statistic), NA_character_))
	
	# check model type & get sample size + case numbers
	model = ""
	if (check_model)  {
		if ("glm" %in% class(x))  if (x$family$family == "gaussian")  {
			model = "Linear model (estimate=coefficient)"
			if (get_r2)  text_out = paste0(text_out, " :: R2=", signif(cor(x$y,predict(x))^2,2))
		}
		if ("glm" %in% class(x))  if (x$family$family == "binomial")  {
			model = "Binomial model (estimate=Odds Ratio)"
			text_out = paste0(text_out, ", Ncases=", as.numeric(table(x$y)[2]))
		}
		if (any(c("coxph") %in% class(x)))  {
			model = "CoxPH model (estimate=Hazard Ratio)"
			text_out = paste0(text_out, ", Nevents=", x$nevent)
		}
		if (any(c("crr","tidycrr") %in% class(x)))  {
			model = "CRR model (estimate=sub-Hazard Ratio)"
			crr_n = table(x$data[,2])
			crr_n2 = crr_n1 = NA
			if (!is.na(crr_n["1"]))  crr_n1 = as.numeric(crr_n["1"])
			if (!is.na(crr_n["2"]))  crr_n2 = as.numeric(crr_n["2"])
			if (is.na(crr_n["2"]))   crr_n2 = 0
			text_out = paste0(text_out, ", Nevents=", crr_n1, ", Ncompeting=", crr_n2)
		}
		text_out = paste0(model, " :: ", text_out)
	}
	
	# exponentiate estimate and CIs?
	if (check_model & !exp & model != "")  exp = TRUE
	if (exp) ret = ret |> dplyr::mutate(estimate=exp(estimate))
	if (exp & (ci | conf.int)) ret = ret |> dplyr::mutate(conf.low=exp(conf.low), conf.high=exp(conf.high))
	
	# print message?
	if (!quiet)  cat(paste0(text_out, "\n"))
	
	# return object
	ret
	
}
