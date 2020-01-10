# Dataset transformations

#' Kendall transformation
#' 
#' @param x Vector or data frame to be Kendall-transformed; allowed feature types are numeric, integer (treated as numeric), ordered factor, logical and unordered factor with two or less levels.
#' \code{NA} and non-finite values are allowed; \code{NaN} is treated as \code{NA}.
#' @return A transformed vector or data frame with transformed columns.
#' @examples
#' kTransform(data.frame(Asc=1:3,Desc=3:1,Vsh=c(2,1,2)))
#' @export
kTransform<-function(x)
 if(is.data.frame(x)) data.frame(.Call(C_kt,x)) else .Call(C_kt,x)
