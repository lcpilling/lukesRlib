#' inv_norm
#'
#' @description Inverse normalize a variable
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
#' inv_norm(x)
#'
#' @export
#'

inv_norm = function(x = stop("Provide a vector")) {
	
	if (!is.numeric(x))  stop("x needs to be numeric")
	
	x_in = qnorm( (rank(x, na.last="keep") - 0.5) / length(na.omit(x)) )
	
	## return vector
	x_in
	
}
