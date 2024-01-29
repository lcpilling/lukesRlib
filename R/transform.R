#' Inverse normalize a variable
#'
#' @description Inverse normalize a variable (force a normal distribution)
#'
#' @return Returns a vector
#'
#' @author Luke Pilling
#'
#' @name inv_norm
#'
#' @param x a numeric vector to be inverse normalized
#'
#' @examples
#' summary(example_data$sbp)
#' sbp_in = inv_norm(example_data$sbp)
#' summary(sbp_in)
#'
#' @export
#'

inv_norm = function(x) {
	if (!is.numeric(x))  stop("x needs to be numeric")
	if (length(x)<500)  warning("x has length<500 -- this may not produce 'nice' output")
	x_in = qnorm( (rank(x, na.last="keep") - 0.5) / length(na.omit(x)) )
	x_in
}

#' Z-transform a variable
#'
#' @description Z-transform a variable. Mean=0, SD=1. Maintains original distribution.
#'
#' @return Returns a vector
#'
#' @author Luke Pilling
#'
#' @name z_trans
#'
#' @param x a numeric vector to be z-transformed
#'
#' @examples
#' summary(example_data$sbp)
#' sbp_in = z_trans(example_data$sbp)
#' summary(sbp_in)
#'
#' @export
#'

z_trans = function(x) {
	if (!is.numeric(x))  stop("x needs to be numeric")
	x_zt = ( ( x - mean(x, na.rm=TRUE) ) / sd(x, na.rm=TRUE) )
	x_zt
}



