#' Get extreme p-value
#'
#' @description Return p-value for any z (or t) statistic as a string
#'
#' @return Returns a string (p-value)
#'
#' @author The internet (anon)
#'
#' @name get_extreme_p
#'
#' @param z a z (or t) statistic
#'
#' @examples
#' z = 50
#' get_extreme_p(z)
#' #>  [1] "2.16e-545"
#'
#' @export
#'

get_extreme_p = function(z) {
	if (!is.numeric(z))  stop("z needs to be numeric")
	log_pvalue   = log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
	log10_pvalue = log_pvalue/log(10) ## from natural log to log10
	mantissa     = 10^(log10_pvalue %% 1)
	exponent     = log10_pvalue %/% 1
	return(sprintf("%1.2fe%d",mantissa,exponent))
}


#' Get -log10 p-value
#'
#' @description Get -log10 p-value from a provided z statistic and n
#'
#' @return Returns a -log 10 p-value
#'
#' @author The internet (anon)
#'
#' @name get_neglog10_p
#'
#' @param z a z (or t) statistic
#' @param n the sample size used to estimate the z (or t) statistic
#'
#' @examples
#' z = 50
#' n = 100000
#' get_neglog10_p(z, n)
#' #>  [1] 537.9851
#'
#' @export
#'

get_neglog10_p = function(z, n) {
	if (!is.numeric(z))  stop("z needs to be numeric")
	if (!is.numeric(n))  stop("n needs to be numeric")
	neglog10_p=-1*(pt(abs(z),df=n,lower.tail=F,log.p=T) + log(2))/log(10)
	neglog10_p
}


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
#'
#' @examples
#' p = 1e-10
#' get_z(p)
#' #>  [1] 6.466951
#'
#' @export
#'
get_z = function(p) {
	if (!is.numeric(p))  stop("p needs to be numeric")
	z = abs(qnorm(p/2,0,1))
	z
}

