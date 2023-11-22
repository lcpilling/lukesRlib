#' Get tidy model output of exposure(s) on outcome(s)
#'
#' @description To easily get tidy model output (from linear or logistic GLM, or CoxPH) for a categorical or continuous exposure, including sample size (and N cases if logistic), outcome, and model info. Idea is to make quick PheWAS trivially easy.
#'
#' For all exposures, it gets the N. For categorical exposures, the N is split by group, and a row is included for the reference category
#'
#' Can provide multiple exposures and/or outcomes 
#'
#' @return Returns a tibble - summary statistics from a model 
#'
#' @author Luke Pilling
#'
#' @name get_assoc
#'
#' @param d A \code{data.frame} or \code{tibble}. The data.
#' @param x A string or vector of strings. The exposure variable name(s), found in `d`. (character)
#' @param y A string. The outcome variable name, found in `d` -- if model is 'coxph' then paste the survival object here e.g., 'Surv(time, event)' where `time` and `event` are variables in `d`. (character)
#' @param z A string. The covariate formula (e.g., " + age + sex"), found in `d`. Default is no covariates. 
#'        \code{default=""} (character)
#' @param model A string. The type of model to perform. Can be "lm", "logistic" or "coxph".
#'        \code{default="lm"} (character)
#' @param af Logical. Is `x` categorical? I.e., include in formula as `as.factor(x)`.
#'        \code{default=FALSE}
#' @param note A string. If you want to include a note like "All", "Males", "C282Y homozygotes" to describe the model or sample.
#'        \code{default=""} (character)
#' @param get_fit Logical. Default is FALSE. Get model fit? (for each model: lm=R2, logistic=McFadden's pseudo-R2, coxph=Harrell's c-statistic).
#'        \code{default=FALSE}
#' @param extreme_ps Logical. Default is FALSE. If p==0 then return "extreme p-values" as strings.
#'        \code{default=FALSE}
#' @param scale_x Logical. Default is FALSE. Apply scale() function to exposure?
#'        \code{default=FALSE}
#' @param scale_y Logical. Default is FALSE. Apply scale() function to outcome?
#'        \code{default=FALSE}
#' @param inv_norm_x Logical. Apply inv_norm() function to exposure?
#'        \code{default=FALSE}
#' @param inv_norm_y Logical. Apply inv_norm() function to outcome?
#'        \code{default=FALSE}
#' @param subset_d Logical. Subset `d` before passing to analysis (much quicker when multiple exposures/outcomes used).
#'        \code{default=TRUE}
#' @param progress Logical. Show progress bar from {purrr} `map()` function (useful when multiple exposures/outcomes provided).
#'        \code{default=TRUE}
#' @param beep Logical. Beep when done,
#'        \code{default=FALSE}
#' @param beep_sound Numeric. Which sound to use? See `beepr` docs,
#'        \code{default=3}
#' @param verbose Logical. Be verbose,
#'        \code{default=FALSE}
#' @param ... Other `tidy_ci()` options
#'
#' @examples
#' # for one outcome, equivalent to `tidy_ci(glm(sbp ~ bmi +age+sex, d=example_data))` - with added `n`
#' get_assoc(x="bmi", y="sbp", z="+age+sex", d=example_data)
#'
#' # categorical exposure, binary outcome, and stratified analysis (with note)
#' #  - note that data can be passed using the pipe if desired
#' example_data |> dplyr::filter(sex==1) |>
#'   get_assoc(x="bmi_cat", y="event", z="+age", model="logistic", af=TRUE, note="Males only")  |> print(width=500)
#'
#' # multiple exposures and/or outcomes - get pseudo R^2
#' x_vars = c("bmi","sbp","dbp","scl")
#' y_vars = c("event","sex")
#' get_assoc(x=x_vars, y=y_vars, z="+age", d=example_data, model="logistic", get_fit=TRUE)  |> print(width=500)
#'
#' @export
#'
get_assoc = function(
	d, 
	x, 
	y, 
	z = "", 
	model = "lm", 
	af = FALSE, 
	note = "", 
	get_fit = FALSE,
	extreme_ps = FALSE,
	scale_x = FALSE, 
	scale_y = FALSE,
	inv_norm_x = FALSE, 
	inv_norm_y = FALSE,
	subset_d = TRUE,
	progress = TRUE,
	beep = FALSE,
	beep_sound = 3,
	verbose = FALSE,
	...
)  {
	
	# check inputs
	if (class(x) != "character")  stop("x needs to be a string or character vector")
	if (class(y) != "character")  stop("y needs to be a string or character vector")
	if (class(z) != "character")  stop("z needs to be a string")
	if (! any(class(d) %in% c("data.frame","tbl","tbl_df")))  stop("d needs to be a data.frame or tibble. Best to explicitly provide inputs using `get_assoc(x=\"BMI\", y=\"diabetes\", z=\"sex\", d=data)` or `data |> get_assoc(x=\"BMI\", y=\"diabetes\", z=\"sex\")`")
	if (! any(model %in% c("lm","logistic","coxph")))         stop("model needs to be 'lm' 'logistic' or 'coxph'")
	
	# verbose?
	if (verbose)  {
		cat("Model = ", model, "\n")
		cat("N exposures = ", length(x), "\n")
		cat("N outcomes = ", length(y), "\n")
		cat("Exposures categorical? ", af, "\n")
		cat("Note: ", note, "\n")
	}
	
	# if doing coxph, extract variable names for checking:
	yy = y
	if (model == "coxph")  yy = stringr::str_replace_all(y, stringr::fixed("Surv("), "") |> 
	                            stringr::str_replace_all(" |\\)", "") |> 
	                            stringr::str_split(",") |> 
	                            purrr::list_c()
	
	# check variables are all in d
	if (any(! x  %in% colnames(d)))  stop("Not all exposure variables are in the provided data")
	if (any(! yy %in% colnames(d)))  stop("Not all outcome variables are in the provided data")
	z_vars = stringr::str_replace_all(z, " |as.factor\\(|\\)", "") |> 
	         stringr::str_split_1(stringr::fixed("+")) |> 
	         purrr::keep(\(x) stringr::str_length(x)>0)
	if (any(! z_vars %in% colnames(d)))  stop("Not all covariate variables are in the provided data")
	if (verbose)  cat("All x, y and z variables are in d\n")
	
	# subset d to just variables used - makes processing much quicker
	if (subset_d)  {
		all_vars = c(x,yy,z_vars)
		d = d[,colnames(d) %in% all_vars]
	}
	if (verbose)  cat("Subsetted data\n")
	
	# check z formula starts with a "+" - if not, add one (unless string is empty)
	if (z != "")  {
		z = stringr::str_replace_all(z, " ", "")
		if (stringr::str_sub(z,start=1,end=1) != "+")  z = paste0("+",z)
	}
	
	# use {purrr} function map2() for analysis -- allows for any combination of exposure/outcome numbers
	if (verbose)  cat("Beginning analysis:\n\n")
	ret = purrr::map2(lukesRlib:::xv(x,y), 
	                  lukesRlib:::yv(x,y), 
	                  \(x,y) lukesRlib:::get_assoc1(x=x, y=y, z=z, d=d, 
	                                                model=model, af=af, note=note, get_fit=get_fit,
	                                                scale_x=scale_x, scale_y=scale_y, inv_norm_x=inv_norm_x, inv_norm_y=inv_norm_y,
	                                                verbose=verbose), 
	                 .progress = progress) |> 
	                 purrr::list_rbind()
	
	# Get extreme p-values if any p-values rounded to zero? 
	if (extreme_ps)  {
		if (verbose)  cat("Getting extreme p-values\n")
		if (any(ret$p.value==0, na.rm=TRUE))  ret = ret |> dplyr::mutate(p.extreme=dplyr::if_else(p.value==0, lukesRlib::get_p_extreme(statistic), NA_character_))
	}
	
	# try to beep?
	if (beep)  try(beepr::beep(beep_sound), silent=TRUE)
	
	# return results
	ret
	
}


#' Internal function to get tidy model output (just exposure of interest, with augmented output)
#' Inherits most options from `get_assoc()`
#' @noRd
get_assoc1 = function(
	x, 
	y, 
	z, 
	d, 
	model = "lm", 
	af = FALSE, 
	note = "", 
	scale_x = FALSE, 
	scale_y = FALSE,
	inv_norm_x = FALSE, 
	inv_norm_y = FALSE,
	get_fit = FALSE,
	verbose = FALSE,
	...
)  {
	
	if (verbose)  cat("Doing analysis of x (", x, ") on y (", y, ")\n")
	
	# if x == y skip this
	if (x != y)  {
		
		# what model are we doing?
		logistic = FALSE
		coxph = FALSE
		if (model == "logistic")  logistic = TRUE
		if (model == "coxph")     coxph = TRUE
		
		# if doing coxph, extract variable names for checking:
		yy = y
		if (coxph)  {
			y12 = stringr::str_replace_all(y, stringr::fixed("Surv("), "") |> stringr::str_replace_all(" |\\)", "") |> stringr::str_split_1(",")
			y1  = y12[1]
			y2  = y12[2]
			
			# if time variable is not numeric try to convert (and give warning)
			if ( ! d |> dplyr::select(!!y1) |> dplyr::pull() |> is.numeric() )  {
				warning("Time variable not numeric -- attempting conversion with `as.numeric()`")
				d = d |> dplyr::mutate(!! rlang::sym(y1) := as.numeric(!! rlang::sym(y1)) )
			}
		} else {
			yy = stringr::str_c("`", yy, "`")  # add backticks to protect variable name in regression formula
		}
		
		# exposure variable - categorical?
		xx = stringr::str_c("`", x, "`")  # add backticks to protect variable name in regression formula
		if (af)  xx = paste0("as.factor(",xx,")")
		
		# scale exposure or outcome?
		if (scale_x & !af)                 xx = paste0("scale(",xx,")")
		if (scale_y & !logistic & !coxph)  yy = paste0("scale(",yy,")")
		
		# inverse normalize exposure or outcome?
		if (inv_norm_x & !af)                d = d |> dplyr::mutate( !! rlang::sym(x) := lukesRlib::inv_norm( !! rlang::sym(x) ) )
		if (inv_norm_y & !logistic & !coxph) d = d |> dplyr::mutate( !! rlang::sym(y) := lukesRlib::inv_norm( !! rlang::sym(y) ) )
		
		# run model
		reg_formula <- paste0(yy, " ~ ", xx, z)
		if (verbose)  cat("Data checked, running model\n - formula:", reg_formula, "\n")
		
		if (!logistic & !coxph)  fit = glm(as.formula(reg_formula), data=d)
		if (logistic)            fit = glm(as.formula(reg_formula), data=d, family=binomial(link="logit"))
		if (coxph)               fit = survival::coxph(as.formula(reg_formula), data=d)
		
		# get tidy output
		res = lukesRlib::tidy_ci(fit, extreme_ps=FALSE, quiet=TRUE, get_r2=FALSE, ...) |> dplyr::filter(grepl(!!x, term))
		
		# include outcome name as first col
		if (coxph)  res = res |> dplyr::mutate(outcome=!!y2) |> dplyr::relocate(outcome)
		if (!coxph) res = res |> dplyr::mutate(outcome=!!y)  |> dplyr::relocate(outcome)
		
		# exposure categorical?
		if (af)  {
			
			x_vals = as.vector(fit$xlevels[[1]])
			
			# categorical exposure - add reference group
			res = rbind(res[1,], res)
			res[1,3:ncol(res)] = NA
			res[1,2] = paste0(x, "-", x_vals[1])
			
			# get sample size - categorical exposure
			if (!coxph)  {
				n = x_vals_n = d |> dplyr::select(!!x, !!y) |> na.omit() |> dplyr::select(!!x) |> table()
				
				# make sure values line up with x labels
				for (ii in 1:length(x_vals))  {
					n[ii] = x_vals_n[x_vals[ii]]
					if (is.na(n[ii]))  n[ii] = 0
				}
				res = res |> dplyr::mutate(n=as.numeric(n))
			}
			
			# if model is logistic, get the cases too
			if (logistic)  {
				n_cases = x_vals_n_cases = d |> dplyr::select(!!x, !!y) |> na.omit() |> dplyr::filter(.data[[y]]==1) |> dplyr::select(!!x) |> table()
				
				# make sure values line up with x labels
				for (ii in 1:length(x_vals))  {
					n_cases[ii] = x_vals_n_cases[x_vals[ii]]
					if (is.na(n_cases[ii]))  n_cases[ii] = 0
				}
				res = res |> dplyr::mutate(n_cases=as.numeric(n_cases))
			}
			
			# if model is coxph, get the cases too
			if (coxph)  {
				n = x_vals_n = d |> dplyr::select(!!x, !!y1, !!y2) |> na.omit() |> dplyr::select(!!x) |> table()
				n_cases = x_vals_n_cases = d |> dplyr::select(!!x, !!y1, !!y2) |> na.omit() |> dplyr::filter(.data[[y2]]==1) |> dplyr::select(!!x) |> table()
				
				# make sure values line up with x labels
				for (ii in 1:length(x_vals))  {
					n[ii] = x_vals_n[x_vals[ii]]
					if (is.na(n[ii]))  n[ii] = 0
					n_cases[ii] = x_vals_n_cases[x_vals[ii]]
					if (is.na(n_cases[ii]))  n_cases[ii] = 0
				}
				res = res |> dplyr::mutate(n=as.numeric(n), n_cases=as.numeric(n_cases))
			}
			
		} else {  # exposure is not categorical
			
			# get sample size - continuous exposure
			if (!coxph)  {
				n = length(fit$y)
				res = res |> dplyr::mutate(n)
			}
			
			# if model is logistic, get the cases too
			if (logistic)  {
				n_cases = length(fit$y[fit$y==1])
				res = res |> dplyr::mutate(n_cases)
			}
			
			# if model is coxph, get the cases too
			if (coxph)  {
				n = fit$n
				n_cases = fit$nevent
				res = res |> dplyr::mutate(n, n_cases)
			}
			
		}
		
		# get model fit?
		if (get_fit)  {
			if (verbose)  cat("Getting fit statistic\n")
			fit_stat = NA
			if (!logistic & !coxph)  {
				# linear model - get R^2
				fit_stat = cor(fit$y,predict(fit))^2
			}
			if (logistic)  {
				# logistic model - get pseudo-R^2 (1 - log(likelihood) / log(likelihood - null model))
				#  - need to make sure null model has same N has full model - exclude any missing X or Z
				z_vars = stringr::str_replace_all(z, " |as.factor\\(|\\)", "") |> 
				         stringr::str_split_1(stringr::fixed("+")) |> 
				         purrr::keep(\(x) stringr::str_length(x)>0)
				if (length(z_vars)==0)  z_vars = yy
				fit_null = glm(paste0(yy, " ~ 1"), data=d |> dplyr::select(yy, x, z_vars) |> na.omit(), family=binomial(link="logit"))
				fit_stat = 1-logLik(fit)/logLik(fit_null)
				fit_stat = fit_stat[1]
			}
			if (coxph)  {
				# CoxPH model - get C-statistic
				fit_stat = concordance(fit)
				fit_stat = fit_stat$concordance
			}
			res = res |> dplyr::mutate(fit_stat)
		}
		
		# modify final bits
		res = res |> dplyr::rename(exposure=term)
		res = res |> dplyr::mutate(model)
		if (nchar(note)>0)  res = res |> dplyr::mutate(note=!!note)
		
		res
	} else {
		if (verbose)  cat("x == y :. skipping (", x, " == ", y, ")\n")
	}
}


#' Get exposure variables for using with `purrr::map2()`
#'
#' @description The {purrr} package provides the wonderful `map()` and `map2()` functions. To get `map2()` to do each comparison need to provide a "long" version of the variable list
#'
#' @return Returns a vector of strings
#'
#' @author Luke Pilling
#'
#' @name xv
#'
#' @param x_vars A vector of strings. The exposure variable names
#' @param y_vars A vector of strings. The outcome variable names
#'
#' @examples
#' x_vars = c("bmi","sbp")
#' y_vars = c("chd","t2d")
#' x_vars2 = xv(x_vars, y_vars)
#' y_vars2 = yv(x_vars, y_vars)
#'
#' @noRd

xv = function(x_vars, y_vars) {
	purrr::map(x_vars, \(x) rep(x, length(y_vars))) |> purrr::list_c()
}


#' Get outcome variables for using with `purrr::map2()`
#'
#' @description The {purrr} package provides the wonderful `map()` and `map2()` functions. To get `map2()` to do each comparison need to provide a "long" version of the variable list
#'
#' @return Returns a vector of strings
#'
#' @author Luke Pilling
#'
#' @name yv
#'
#' @param x_vars A vector of strings. The exposure variable names
#' @param y_vars A vector of strings. The outcome variable names
#'
#' @examples
#' x_vars = c("bmi","sbp")
#' y_vars = c("chd","t2d")
#' x_vars2 = xv(x_vars, y_vars)
#' y_vars2 = yv(x_vars, y_vars)
#'
#' @noRd

yv = function(x_vars, y_vars) {
	purrr::map(x_vars, \(x) c(y_vars)) |> purrr::list_c()
}


