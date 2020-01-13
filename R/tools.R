#' Join factors
#'
#' Convenience function for joining factors.
#' 
#' @param ... One or more features to merge. 
#' Given as single vectors or data.frames. 
#' Accepted feature types are factor (preferred), booleans, integers (treated as categorical) or reals (which undergo automatic categorisation). 
#' \code{NA}s are not allowed.
#' @return Joint factor, with levels \code{l1} to \code{l<n>}.
#' Vacant combinations are dropped. 
#' @note You can pass a single vector to this function to see how praznik interprets it.
#' @examples
#' joinf(c(T,F,T,F),c(T,T,F,F))
#' @export
joinf<-function(...)
 .Call(C_join,data.frame(...))
