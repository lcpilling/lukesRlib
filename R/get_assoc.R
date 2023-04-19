#' Get tidy model output of exposure(s) on outcome(s)
#'
#' @description To easily get tidy model output for a categorical or continuous exposure, including sample size (and N cases if logistic), outcome, and model info. Idea is to make quick loops easy.
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
#' @param x A string or vector of strings. The exposure variable name(s), found in `d`
#' @param y A string or vector of strings. The outcome variable name(s), found in `d`
#' @param z A string. The covariate formula (e.g., " + age + sex"), found in `d`
#' @param d A data.frame or tibble. The data
#' @param subset_d Logical. Default is TRUE. Subset `d` before passing to analysis (much quicker when multiple exposures/outcomes used)
#' @param progress Logical. Default is TRUE. Show progress bar from {purrr} `map()` function (useful when multiple exposures/outcomes provided)
#' @param logistic Logical. Default is FALSE. Include `family=binomial(link="logit")` in `glm()`?
#' @param af Logical. Default is FALSE. Is `x` categorical? I.e., include in formula as `as.factor(x)`
#' @param note A string. If you want to include a note like "All", "Males", "C282Y homozygotes" to describe the model or sample.
#' @param scale_x Logical. Default is FALSE. Apple scale() function to exposure?
#' @param scale_y Logical. Default is FALSE. Apple scale() function to outcome?
#' @param extreme_ps Logical. Default is TRUE. If p==0 then return "extreme p-values" as strings.
#' @param ... Other `tidy_ci()` options
#'
#' @examples
#' # for one outcome, equivalent to `tidy_ci(glm(weight ~ height +age+sex, d=ukb))` - with added `n`
#' get_assoc(y="weight", x="height", z="+age+sex", d=ukb)
#'
#' # categorical exposure, binary outcome, and stratified analysis (with note)
#' get_assoc(y="chd", x="smoking_status", z="+age", d=ukb |> filter(sex==1), logistic=TRUE, af=TRUE, note="Males only")
#'
#' # multiple exposures and/or outcomes
#' x_vars = c("bmi","ldl","sbp_0_avg")
#' y_vars = c("chd","t2d")
#' get_assoc(x=x_vars, y=y_vars, z="+age+sex", d=ukb, logistic=TRUE)
#'
#' @export
#'
get_assoc = function(x, y, z, d, 
                     subset_d = TRUE,
                     progress = TRUE,
                     logistic = FALSE, 
                     af = FALSE, 
                     note = "", 
                     scale_x = FALSE, 
                     scale_y = FALSE,
                     extreme_ps = TRUE,
                     ...)  {
	
	# check inputs
	if (class(x) != "character")  stop("x needs to be a string or character vector")
	if (class(y) != "character")  stop("y needs to be a string or character vector")
	if (class(z) != "character")  stop("z needs to be a string")
	if (! any(class(d) %in% c("data.frame","tbl","tbl_df")))  stop("d needs to be a data.frame or tibble")
	
	# check variables are all in d
	if (any(! x %in% colnames(d)))  stop("Not all exposure variables are in the provided data")
	if (any(! y %in% colnames(d)))  stop("Not all exposure variables are in the provided data")
	z_vars = stringr::str_replace_all(z, " |as.factor\\(|\\)", "") |> 
	         stringr::str_split_1(stringr::fixed("+")) |> 
	         purrr::keep(\(x) stringr::str_length(x)>0)
	if (any(! z_vars %in% colnames(d)))  stop("Not all covariate variables are in the provided data")
	
	# subset d to just variables used - makes processing much quicker
	if (subset_d)  {
		all_vars = c(x,y,z_vars)
		d = d[,colnames(d) %in% all_vars]
	}
	
	# use {purrr} function map2() for analysis -- allows for any combination of exposure/outcome numbers
	ret = purrr::map2(lukesRlib::xv(x,y), lukesRlib::yv(x,y), 
	                  \(x,y) lukesRlib::get_assoc1(x=x, y=y, z=z, d=d, logistic=logistic, af=af, note=note, scale_x=scale_x, scale_y=scale_y), 
	                 .progress = progress) |> purrr::list_rbind()
	
	# Get extreme p-values if any p-values rounded to zero? 
	if (extreme_ps)  if (any(ret$p.value==0, na.rm=TRUE))  ret = ret |> dplyr::mutate(p.extreme=dplyr::if_else(p.value==0, lukesRlib::get_p_extreme(statistic), NA_character_))
	
	# return results
	ret
	
}



#' Get tidy model output (just exposure of interest, with augmented output)
#'
#' @description To easily get tidy model output for a categorical or continuous exposure, including sample size (and N cases if logistic), outcome, and model info. Idea is to make quick loops easy.
#'
#' For all exposures, it gets the N. For categorical exposures, the N is split by group, and a row is included for the reference category
#'
#' !!!! Better to just use `get_assoc()` though, which calls this function but after checking variables and subsetting etc.
#'
#' @return Returns a tibble - summary statistics from a model 
#'
#' @author Luke Pilling
#'
#' @name get_assoc1
#'
#' @param x A string. The exposure variable name, found in `d`
#' @param y A string. The outcome variable name, found in `d`
#' @param z A string. The covariate formula (e.g., " + age + sex"), found in `d`
#' @param d A data.frame or tibble. The data
#' @param logistic Logical. Default is FALSE. Include `family=binomial(link="logit")` in `glm()`?
#' @param af Logical. Default is FALSE. Is `x` categorical? I.e., include in formula as `as.factor(x)`
#' @param note A string. If you want to include a note like "All", "Males", "C282Y homozygotes" to describe the model or sample.
#' @param scale_x Logical. Default is FALSE. Apple scale() function to exposure?
#' @param scale_y Logical. Default is FALSE. Apple scale() function to outcome?
#' @param ... Other `tidy_ci()` options
#'
#' @examples
#' # for one outcome, equivalent to `tidy_ci(glm(weight ~ height +age+sex, d=ukb))` - with added `n`
#' get_assoc1(y="weight", x="height", z="+age+sex", d=ukb)
#'
#' # categorical exposure, binary outcome, and stratified analysis (with note)
#' get_assoc1(y="chd", x="smoking_status", z="+age", d=ukb |> filter(sex==1), logistic=TRUE, af=TRUE, note="Males only")
#'
#' # multiple exposures on single outcome, then combine output
#' x_vars = c("bmi","ldl","sbp_0_avg")
#' res = do.call(rbind, lapply(x_vars, get_assoc1, y="chd", z="+age+sex", d=ukb, logistic=TRUE))
#'
#' # one exposure on multiple outcomes, then combine output
#' y_vars = c("bmi","ldl","sbp_0_avg")
#' res = do.call(rbind, lapply(y_vars, function(y) get_assoc1(x="smoking_status", y=y, z="+age+sex", d=ukb, af=TRUE)))
#'
#' @export
#'

get_assoc1 = function(x, y, z, d, 
                     logistic = FALSE, 
                     af = FALSE, 
                     note = "", 
                     scale_x = FALSE, 
                     scale_y = FALSE,
                     ...)  {

	# check inputs
	if (class(x) != "character")  stop("x needs to be a string, a variable name in d")
	if (class(y) != "character")  stop("y needs to be a string, a variable name in d")
	if (class(z) != "character")  stop("z needs to be a string of variable names in d")
	if (! any(class(d) %in% c("data.frame","tbl","tbl_df")))  stop("d needs to be a data.frame or tibble")

	if (! x %in% colnames(d))  stop("x needs to be a variable in d")
	if (! y %in% colnames(d))  stop("y needs to be a variable in d")

	# exposure variable - categorical?
	xx = x
	if (af)  xx = paste0("as.factor(",x,")")
	
	# scale exposure or outcome?
	yy = y
	if (scale_x & !af)  xx = paste0("scale(",x,")")
	if (scale_y & !logistic)  yy = paste0("scale(",y,")")

	# run model
	if (!logistic)  {
		model = "lm"
		fit = glm(paste0(yy, " ~ ", xx, z), data=d)
	} else {
		model = "logistic"
		fit = glm(paste0(yy, " ~ ", xx, z), data=d, family=binomial(link="logit"))
	}

	# get tidy output - include outcome name as first col
	res = lukesRlib::tidy_ci(fit, extreme_ps=FALSE, quiet=TRUE, get_r2=FALSE, ...) |> dplyr::filter(grepl(!!x, term))
	res = res |> dplyr::mutate(outcome=!!y) |> dplyr::relocate(outcome)

	# exposure categorical?
	if (af)  {

		x_vals = as.vector(fit$xlevels[[1]])

		# categorical exposure - add reference group
		res = rbind(res[1,], res)
		res[1,3:ncol(res)] = NA
		res[1,2] = paste0(x, "-", x_vals[1])

		# get sample size - categorical exposure
		n = x_vals_n = table(d |> dplyr::select(!!x, !!y) |> na.omit() |> dplyr::select(!!x))
		for (ii in 1:length(x_vals))  {
			n[ii] = x_vals_n[x_vals[ii]]
			if (is.na(n[ii]))  n[ii] = 0
		}
		res = res |> dplyr::mutate(n=as.numeric(n))

		# if model is logistic, get the cases too
		if (model == "logistic")  {
			n_cases = x_vals_n_cases = table(d |> dplyr::select(!!x, !!y) |> na.omit() |> dplyr::filter(.data[[y]]==1) |> dplyr::select(!!x))
			for (ii in 1:length(x_vals))  n_cases[ii] = x_vals_n_cases[x_vals[ii]]
			res = res |> dplyr::mutate(n_cases=as.numeric(n_cases))
		}

	} else {  # exposure is not categorical

		# get sample size - continuous exposure
		n = length(fit$y)
		res = res |> dplyr::mutate(n)
		
		# if model is logistic, get the cases too
		if (model == "logistic")  {
			n_cases = length(fit$y[fit$y==1])
			res = res |> dplyr::mutate(n_cases)
		}

	}

	# modify final bits
	res = res |> dplyr::rename(exposure=term)
	res = res |> dplyr::mutate(model)
	if (nchar(note)>0)  res = res |> dplyr::mutate(note=!!note)

	res
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
#' @export
#'

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
#' @export
#'

yv = function(x_vars, y_vars) {
	purrr::map(x_vars, \(x) c(y_vars)) |> purrr::list_c()
}


