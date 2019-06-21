# Dataset transformations

#TODO: Rewrite this in C

#' Kendall transformation
#' 
#' @param x Vector or data frame to be Kendall-transformed.
#' @export
kTransform<-function(x){
 if(is.data.frame(x))
  return(data.frame(lapply(x,kTransform)))
 
 #Prepare pairs
 mask<-expand.grid(a=1:length(x),b=1:length(x))
 mask[mask$a<mask$b,]->mask
 xa<-x[mask$a]
 xb<-x[mask$b]
 
 if(is.numeric(x)||is.ordered(x))
  return(factor(ifelse(xa>xb,'>',ifelse(xa<xb,'<','='))))
 
 #TODO: This requires research
 factor(ifelse(xa==xb,'=','<>'))
 # or
 #factor(ifelse(xa==xb,sprintf("=%s",xa),'<>'))
 # or
 # factor(sprintf("%s_%s",xa,xb)) #TODO: buggy

}
