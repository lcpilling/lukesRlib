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
#' @param af Logical. Is `x` categorical? I.e., include in formula like `haven::as_factor(x)`.
#'        To use base R `as.factor()` instead set option "af_base" to TRUE
#'        \code{default=FALSE}
#' @param af_base Logical. Use base R `as.factor()` instead of `haven::as_factor(x)`?
#'        \code{default=FALSE}
#' @param note A string. If you want to include a note like "All", "Males", "C282Y homozygotes" to describe the model or sample.
#'        \code{default=""} (character)
#' @param get_fit Logical. Default is FALSE. Get model fit? (for each model: lm=R2, logistic=McFadden's pseudo-R2, coxph=Harrell's c-statistic).
#'        \code{default=FALSE}
#' @param extreme_ps Logical. Default is FALSE. If p==0 then return "extreme p-values" as strings.
#'        \code{default=FALSE}
#' @param include_formula Logical. Default is FALSE. Include the regression formula in the output?
#'        \code{default=FALSE}
#' @param scale_x Logical. Default is FALSE. Apply scale() function to exposure?
#'        \code{default=FALSE}
#' @param scale_y Logical. Default is FALSE. Apply scale() function to outcome?
#'        \code{default=FALSE}
#' @param inv_norm_x Logical. Apply inv_norm() function to exposure?
#'        \code{default=FALSE}
#' @param inv_norm_y Logical. Apply inv_norm() function to outcome?
#'        \code{default=FALSE}
#' @param winsorize_x Logical. Apply Winzorization to exposure?
#'        \code{default=FALSE}
#' @param winsorize_y Logical. Apply Winzorization to outcome?
#'        \code{default=FALSE}
#' @param winsorize_n Numeric. Standard deviations from the mean to Winzorize. I.e., participants with values beyond this N will be set to N.
#'        \code{default=5}
#' @param return_all_terms Logical. Return estimates for all independent variables (terms) in the model? If TRUE, includes an adddition column 'term'
#'        \code{default=FALSE}
#' @param interacts_with A string. A variable found in `d`. Will add to regression formula like `x*i` and catch output
#'        \code{default=""} (character)
#' @param progress Logical. Show progress bar from {purrr} `map()` function (useful when multiple exposures/outcomes provided).
#'        \code{default=TRUE}
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
#' # if desired, can also return estimates for other independent variables (terms) in the model
#' x_vars = c("bmi","sbp","dbp")
#' get_assoc(x=x_vars, y="event", z="+age+sex", d=example_data, model="logistic", return_all_terms=TRUE)  |> print(width=500)
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
	af_base = FALSE, 
	note = "", 
	get_fit = FALSE,
	extreme_ps = FALSE,
	include_formula = FALSE,
	scale_x = FALSE, 
	scale_y = FALSE,
	inv_norm_x = FALSE, 
	inv_norm_y = FALSE,
	winsorize_x = FALSE, 
	winsorize_y = FALSE,
	winsorize_n = 5,
	return_all_terms = FALSE,
	interacts_with = "",
	progress = TRUE,
	verbose = FALSE,
	...
)  {
	
	v <- packageVersion("lukesRlib")
	cli::cli_alert_info("lukesRlib v{v}")
	start_time <- Sys.time()
	
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
		
		if (scale_x) cat("Scaling X values\n")
		if (scale_y) cat("Scaling Y values\n")
		if (inv_norm_x) cat("Rank-based inverse normal transforming X values\n")
		if (inv_norm_y) cat("Rank-based inverse normal transforming Y values\n")
		if (winsorize_x) cat("Winsorizing X values (values >", winsorize_n, " SDs from mean are truncated)\n")
		if (winsorize_y) cat("Winsorizing Y values (values >", winsorize_n, " SDs from mean are truncated)\n")
	}
	
	# if doing coxph, extract variable names for checking:
	yy = y
	if (model == "coxph")  yy = stringr::str_replace_all(y, stringr::fixed("Surv("), "") |> 
	                            stringr::str_replace_all(" |\\)", "") |> 
	                            stringr::str_split(",") |> 
	                            purrr::list_c()
	
	# check variables are all in d
	if (any(! x  %in% colnames(d)))  {
		cat("!! Not all exposure variables are in the provided data\n")
		missing <- x[!x %in% colnames(d)]
		if (!verbose & length(missing>6))  {
			cat("Here are the first 6 that are missing...\n")
			return(missing[1:6])
		}
		if (verbose)  print(missing)
		stop()
	}
	if (any(! yy %in% colnames(d)))  {
		cat("!! Not all outcome variables are in the provided data\n")
		missing <- yy[!yy %in% colnames(d)]
		if (!verbose & length(missing>6))  {
			cat("Here are the first 6 that are missing...\n")
			return(missing[1:6])
		}
		if (verbose)  return(missing)
		stop()
	}
	z_vars = stringr::str_replace_all(z, " |as.factor\\(|haven::|\\)", "") |> 
	         stringr::str_split_1("\\+|\\*") |> 
	         purrr::keep(\(x) stringr::str_length(x)>0)
	if (any(! z_vars %in% colnames(d)))  {
		cat("!! Not all covariate variables are in the provided data\n")
		missing <- z_vars[!z_vars %in% colnames(d)]
		if (!verbose & length(missing>6))  {
			cat("Here are the first 6 that are missing...\n")
			return(missing[1:6])
		}
		if (verbose)  print(missing)
		stop()
	}
	if (verbose)  cat("All x, y and z variables are in d\n")
	
	# doing interactions?
	if (stringr::str_length(interacts_with)>0)  {
		if (any(! interacts_with %in% colnames(d)))  stop("Interaction variable is not in the provided data")
		if (verbose)  cat("Getting interactions with `", interacts_with, "`\n")
	}
	
	# subset d to just variables used - makes processing much quicker
	all_vars = c(x,yy,z_vars)
	if (stringr::str_length(interacts_with)>0)  all_vars = c(all_vars, interacts_with)
	d = d |> dplyr::select(dplyr::all_of(all_vars))
	if (verbose)  cat("Subsetted data\n")
	
	# check z formula starts with a "+" - if not, add one (unless string is empty)
	if (z != "")  {
		z = stringr::str_replace_all(z, " ", "")
		if (stringr::str_sub(z,start=1,end=1) != "+")  z = paste0("+",z)
	}
	
	# use {purrr} function map2() for analysis -- allows for any combination of exposure/outcome numbers
	if (verbose)  cat("Beginning analysis:\n")
	
	# get x and y var lists
	xv_vars <- lukesRlib:::xv(x,y)
	yv_vars <- lukesRlib:::yv(x,y)
	xv_vars_n <- length(xv_vars)
	cli::cli_alert("Getting {xv_vars_n} association{?s}")
	ret = purrr::map2(xv_vars, 
	                  yv_vars, 
	                  \(x,y) lukesRlib:::get_assoc1(x=x, y=y, z=z, d=d, 
	                                                model=model, af=af, af_base=af_base, note=note, get_fit=get_fit, include_formula=include_formula,
	                                                scale_x=scale_x, scale_y=scale_y, 
	                                                inv_norm_x=inv_norm_x, inv_norm_y=inv_norm_y, 
	                                                winsorize_x=winsorize_x, winsorize_y=winsorize_y, winsorize_n=winsorize_n,
	                                                return_all_terms=return_all_terms, interacts_with=interacts_with,
	                                                verbose=verbose), 
	                 .progress = progress) |> 
	                 purrr::list_rbind()
	
	# Get extreme p-values if any p-values rounded to zero? 
	if (extreme_ps)  {
		if (verbose)  cat("\nGetting extreme p-values\n")
		if (any(ret$p.value==0, na.rm=TRUE))  ret = ret |> dplyr::mutate(p.extreme=dplyr::if_else(p.value==0, lukesRlib::get_p_extreme(statistic), NA_character_))
	}
	
	# finished!
	cli::cli_alert_success(c("Finished. Time taken: ", "{prettyunits::pretty_sec(as.numeric(difftime(Sys.time(), start_time, units=\"secs\")))}."))
	
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
	af_base = FALSE, 
	note = "", 
	get_fit = FALSE,
	include_formula=FALSE,
	scale_x = FALSE, 
	scale_y = FALSE,
	inv_norm_x = FALSE, 
	inv_norm_y = FALSE,
	winsorize_x = FALSE, 
	winsorize_y = FALSE,
	winsorize_n = 5,
	return_all_terms = FALSE,
	interacts_with = "",
	verbose = FALSE,
	...
)  {
	
	if (verbose)  cat("\nDoing analysis of x (", x, ") on y (", y, ")\n")
	
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
			if ( ! d |> dplyr::select(dplyr::all_of(y1)) |> dplyr::pull() |> is.numeric() )  {
				warning("Time variable not numeric -- attempting conversion with `as.numeric()`")
				d = d |> dplyr::mutate(!! rlang::sym(y1) := as.numeric(!! rlang::sym(y1)) )
			}
		} else {
			yy = stringr::str_c("`", yy, "`")  # add backticks to protect variable name in regression formula
		}
		
		# add backticks to protect variable name in regression formula
		xx = stringr::str_c("`", x, "`")  

		# start formula
		this_formula_y <- yy
		this_formula_x <- xx
		
		# outcome: scale or inverse normalize?
		if (scale_y & !logistic & !coxph)  {
			yy <- paste0("scale(",yy,")")
			this_formula_y <- yy
		}
		if (inv_norm_y & !logistic & !coxph)  {
			d <- d |> dplyr::mutate( !! rlang::sym(y) := lukesRlib::inv_norm( !! rlang::sym(y) ) )
			this_formula_y <- stringr::str_c("lukesRlib::inv_norm(", yy, ")")
		}
		
		# exposure: scale or inverse normalize?
		if (scale_x & !af)  {
			xx = paste0("scale(",xx,")")
			this_formula_x <- xx
		}
		if (inv_norm_x & !af)  {
			d = d |> dplyr::mutate( !! rlang::sym(x) := lukesRlib::inv_norm( !! rlang::sym(x) ) )
			this_formula_x <- stringr::str_c("lukesRlib::inv_norm(", xx, ")")
		}
		
		# winsorizing exposure or outcome?
		if (winsorize_x & !af)  {
			v = d |> dplyr::select( !! rlang::sym(x) ) |> dplyr::pull()
			min = mean(v, na.rm=T) - winsorize_n*sd(v, na.rm=T)
			max = mean(v, na.rm=T) + winsorize_n*sd(v, na.rm=T)
			d = d |> dplyr::mutate( !! rlang::sym(x) := dplyr::case_when( 
				!! rlang::sym(x) > max ~ max,
				!! rlang::sym(x) < min ~ min,
				TRUE ~ !! rlang::sym(x)
			) )
		}
		if (winsorize_y & !logistic & !coxph)  {
			v = d |> dplyr::select( !! rlang::sym(y) ) |> dplyr::pull()
			min = mean(v, na.rm=T) - winsorize_n*sd(v, na.rm=T)
			max = mean(v, na.rm=T) + winsorize_n*sd(v, na.rm=T)
			d = d |> dplyr::mutate( !! rlang::sym(y) := dplyr::case_when( 
				!! rlang::sym(y) > max ~ max,
				!! rlang::sym(y) < min ~ min,
				TRUE ~ !! rlang::sym(y)
			) )
		}
		
		# exposure variable - categorical?
		#   if using haven, mutate the actual data
		if (af & af_base)  {
			xx = stringr::str_c("as.factor(",xx,")")
			this_formula_x <- xx
		}
		if (af & !af_base)  {
			d = d |> dplyr::mutate( !! rlang::sym(x) := droplevels( haven::as_factor( !! rlang::sym(x) ) ) )
			xx = stringr::str_c("haven::as_factor(",xx,")")
			this_formula_x <- xx
		}
		
		# put together actual regression formula
		reg_formula <- stringr::str_c(yy, " ~ ", xx)
		if (stringr::str_length(interacts_with)>0)  reg_formula <- stringr::str_c(reg_formula, "*", interacts_with)
		reg_formula <- stringr::str_c(reg_formula, z)
		if (coxph)  reg_formula <- stringr::str_c("survival::", reg_formula)
		
		# put together dummy formula if required for output
		if (include_formula)  {
			this_formula <- stringr::str_c(this_formula_y, " ~ ", this_formula_x)
			if (coxph)  this_formula <- stringr::str_c("survival::", this_formula)
			if (stringr::str_length(interacts_with)>0)  this_formula <- stringr::str_c(this_formula, "*", interacts_with)
			this_formula <- stringr::str_c(this_formula, z)
			this_formula <- stringr::str_replace_all(this_formula, stringr::fixed("+"), " + ")
		}

		if (verbose)  {
			cat("Data checked, running model\n")
			cat("- real formula:", reg_formula, "\n")
			cat("- dummy formula:", this_formula, "\n")
		}
		
		# run model
		if (!logistic & !coxph)  fit = glm(as.formula(reg_formula), data=d)
		if (logistic)            fit = glm(as.formula(reg_formula), data=d, family=binomial(link="logit"))
		if (coxph)               fit = survival::coxph(as.formula(reg_formula), data=d)
		
		if (verbose)  cat("Model finished - processing\n")
		
		# get tidy output
		res_all = lukesRlib::tidy_ci(fit, extreme_ps=FALSE, quiet=TRUE, get_r2=FALSE, ...)
		if (verbose)  cat("Got tidy output\n")
		
		# Filter to just the exposure variable 
		if (x %in% res_all$term)  {
			res = dplyr::filter(res_all, term == !!x)
		}  else  {
			# need to make sure not capturing covariates with similar names (e.g., if x="age" do not get covariate "percentage")
			# if no exact matches then it is categorical or scaled. Grep for x"-"
			res = dplyr::filter(res_all, stringr::str_detect(term, stringr::str_c(!!x, "-")))
			# exclude if contains ":" as this is an interaction term and will be caught later
			res = dplyr::filter(res, stringr::str_detect(term, stringr::fixed(":"), negate=TRUE))
		}
		
		# include outcome name as first col
		if (coxph)  res = res |> dplyr::mutate(outcome=!!y2) |> dplyr::relocate(outcome)
		if (!coxph) res = res |> dplyr::mutate(outcome=!!y)  |> dplyr::relocate(outcome)
		
		# exposure categorical?
		if (af)  {
			
			if (verbose)  cat("Categorical exposure\n")
			
			x_vals = as.vector(fit$xlevels[[1]])
			if (verbose)  {
				cat("x_vals\n")
				print(x_vals)
			}
			
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
			
			if (verbose)  cat("Continuous exposure\n")
			
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
			fit_stat_se = NA
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
				fit_null = glm(paste0(yy, " ~ 1"), 
				               data=d |> dplyr::select(dplyr::all_of(c(stringr::str_replace_all(yy, "`", ""), x, z_vars))) |> na.omit(), 
				               family=binomial(link="logit"))
				fit_stat = 1-logLik(fit)/logLik(fit_null)
				fit_stat = fit_stat[1]
			}
			if (coxph)  {
				# CoxPH model - get C-statistic
				fit_stat = concordance(fit)
				fit_stat_se = sqrt(fit_stat$var)
				fit_stat = fit_stat$concordance
			}
			res = res |> dplyr::mutate(fit_stat, fit_stat_se)
		}
		
		# rename 'term' to 'exposure'
		res = res |> dplyr::rename(exposure=term)
		
		# include other independent variables in output?
		if (return_all_terms)  {
			
			# add 'term' to `res` - replace var in 'exposure' with just variable name
			res = res |> dplyr::mutate(term=exposure, exposure=!!x) |> dplyr::relocate(term, .after=exposure)
			
			# remove exposure from `res_all`
			res_all = dplyr::filter(res_all, ! term %in% res$term)
			
			# make sure includes other relevant cols (blank, as required)
			res_all = res_all |> 
				dplyr::mutate(exposure=!!x) |> 
				dplyr::relocate(exposure, .before=term) |>
				dplyr::mutate(outcome=unique(res$outcome)) |>
				dplyr::relocate(outcome) |>
				dplyr::mutate(n=NA)
			if ("n_cases" %in% colnames(res))  res_all = dplyr::mutate(res_all, n_cases=NA)
			if (get_fit)  res_all = dplyr::mutate(res_all, fit_stat=NA, fit_stat_se=NA)
			
			# append to `res`
			res = rbind(res, res_all)
			
		}
		
		# include interaction term as output? Not required if including all output
		if (!return_all_terms & stringr::str_length(interacts_with)>0)  {
			
			# add 'term' to `res` - replace var in 'exposure' with just variable name
			res = res |> dplyr::mutate(term=exposure, exposure=!!x) |> dplyr::relocate(term, .after=exposure)
			
			# identify interaction rows
			res_sub = dplyr::filter(res_all, stringr::str_detect(term, !!x) & stringr::str_detect(term, !!interacts_with))
			
			# make sure includes other relevant cols (blank, as required)
			res_sub = res_sub |> 
				dplyr::mutate(exposure=!!x) |> 
				dplyr::relocate(exposure, .before=term) |>
				dplyr::mutate(outcome=unique(res$outcome)) |>
				dplyr::relocate(outcome) |>
				dplyr::mutate(n=NA)
			if ("n_cases" %in% colnames(res))  res_sub = dplyr::mutate(res_sub, n_cases=NA)
			if (get_fit)  res_sub = dplyr::mutate(res_sub, fit_stat=NA, fit_stat_se=NA)
			
			# append to `res`
			res = rbind(res, res_sub)
			
		}
		
		# add 'model', 'note' columns and 'formula' cols as required
		res = res |> dplyr::mutate(model)
		if (nchar(note)>0)    res = res |> dplyr::mutate(note=!!note)
		if (include_formula)  res = res |> dplyr::mutate(formula=!!this_formula)
		
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


