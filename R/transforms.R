# Dataset transformations

#' Kendall transformation
#' 
#' @param x Vector or data frame to be Kendall-transformed.
#' @examples
#' kTransform(data.frame(Asc=1:3,Desc=3:1,Vsh=c(2,1,2)))
#' @export
kTransform<-function(x)
 if(is.data.frame(x)) data.frame(.Call(C_kt,x)) else .Call(C_kt,x)
