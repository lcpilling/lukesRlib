#' Get p-value from Z-statistic
#'
#' @description Return p-value for a z (or t) statistic
#'
#' @return Returns a p-value (numeric)
#'
#' @author The internet (anon)
#'
#' @name get_p
#'
#' @param z a z (or t) statistic
#'
#' @examples
#' z = 10
#' get_p(z)
#' #>  [1] 1.523971e-23
#'
#' @export
#'
get_p = function(z) {
	if (!is.numeric(z))  stop("z needs to be numeric")
	p = 2*pnorm(-abs(z))
	p
}


#' Get Z-statistic from p-value
#'
#' @description Return a z-statistic from a given p-value
#'
#' @return Returns a Z-statistic (numeric)
#'
#' @author The internet (anon)
#'
#' @name get_z
#'
#' @param p a p-value
#' @param two_sided Logical, default=TRUE. Return two-sided z-value?
#'
#' @examples
#' p = 1e-10
#' get_z(p)
#' #>  [1] 6.466951
#'
#' @export
#'
get_z = function(p, two_sided=TRUE) {
	if (!is.numeric(p))  stop("p needs to be numeric")
	if (two_sided)  {
		z = abs(qnorm(p/2))
	}  else  {
		z = abs(qnorm(p))
	}
	z
}


#' Get extreme p-value
#'
#' @description Return p-value for any z (or t) statistic as a string
#'
#' @return Returns a string (p-value)
#'
#' @author The internet (anon)
#'
#' @name get_p_extreme
#'
#' @param z a z (or t) statistic
#'
#' @examples
#' z = 50
#' get_p_extreme(z)
#' #>  [1] "2.16e-545"
#'
#' @export
#'

get_p_extreme = function(z) {
	if (!is.numeric(z))  stop("z needs to be numeric")
	if (!is.na(z)) {
		log_pvalue   = log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
		log10_pvalue = log_pvalue/log(10) ## from natural log to log10
		mantissa     = 10^(log10_pvalue %% 1)
		exponent     = log10_pvalue %/% 1
		return(sprintf("%1.2fe%d",mantissa,exponent))
	} else {
		return(NA)
	}
}


#' Get -log10 p-value (using Z only)
#'
#' @description Get -log10 p-value from a provided z statistic
#'
#' @return Returns a -log 10 p-value
#'
#' @author The internet (anon)
#'
#' @name get_p_neglog10
#'
#' @param z a z (or t) statistic
#'
#' @examples
#' z = 50
#' get_p_neglog10(z)
#' #>  [1] 544.0977
#'
#' @export
#'

get_p_neglog10 = function(z) {
	if (!is.numeric(z))  stop("z needs to be numeric")
	neglog10_p=-( 2 + pnorm(-abs( z ), log.p=T) ) / log(10)
	neglog10_p
}


#' Get -log10 p-value (using Z & N)
#'
#' @description Get -log10 p-value from a provided z statistic and n
#'
#' @return Returns a -log 10 p-value
#'
#' @author The internet (anon)
#'
#' @name get_p_neglog10_n
#'
#' @param z a z (or t) statistic
#' @param n the sample size used to estimate the z (or t) statistic
#'
#' @examples
#' z = 50
#' n = 100000
#' get_p_neglog10_n(z, n)
#' #>  [1] 537.9851
#'
#' @export
#'

get_p_neglog10_n = function(z, n) {
	if (!is.numeric(z))  stop("z needs to be numeric")
	if (!is.numeric(n))  stop("n needs to be numeric")
	neglog10_p=-1*(pt(abs(z),df=n,lower.tail=F,log.p=T) + log(2))/log(10)
	neglog10_p
}


#' Get Standard Error from Confidence Intervals
#'
#' @description Return SE from CIs
#'
#' @return Returns SE (numeric)
#'
#' @author The internet (anon)
#'
#' @name get_se
#'
#' @param lci the lower confidence interval
#' @param uci the lower confidence interval
#' @param log log the CI values? (default=FALSE)
#' @param denominator the standard error of the sample mean (default=3.92... well actually, 3.919928 from `2*get_z(0.05)`)
#'
#' @examples
#' lci = 0.1
#' uci = 0.3
#' get_se(lci, uci)
#' #>  [1] 0.05102041
#'
#' @export
#'
get_se = function(lci, uci, log=FALSE, denominator=3.919928) {
	if (!is.numeric(lci))  stop("`lci` needs to be numeric")
	if (!is.numeric(uci))  stop("`uci` needs to be numeric")
	if (lci > uci)         stop("`uci` needs to be greater than `lci`")
	if (log)  {
		lci = log(lci)
		uci = log(uci)
	}
	se = (uci - lci) / denominator
	se
}
