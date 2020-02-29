# Feature scorers

#' Mutual information scores
#'
#' Calculates mutual information between each feature and the decision, that is
#' \deqn{I(X,Y).}
#' @template input
#' @template y
#' @return A numerical vector with mutual information scores, with names copied from \code{X}.
#' @examples
#' miScores(iris[,-5],iris$Species)
#' @export
miScores<-function(X,Y,threads=0)
 .Call(C_mi,X,Y,as.integer(threads))

#' Mutual information matrix
#'
#' Calculates mutual information between each two features, that is
#' \deqn{I(X_i,X_j).}
#' @template input
#' @template matrix
#' @examples
#' miMatrix(iris)
#' @export
miMatrix<-function(X,zeroDiag=TRUE,threads=0)
 .Call(C_miMatrix,X,as.logical(zeroDiag),as.integer(threads))

#' Normalised mutual information matrix
#'
#' Calculates normalised mutual information between each two features, that is
#' \deqn{\frac{I(X_i,X_j)}{H(X_i,X_j)}.}
#' @template input
#' @template matrix
#' @examples
#' nmiMatrix(iris)
#' @export
nmiMatrix<-function(X,zeroDiag=TRUE,threads=0)
 .Call(C_nmiMatrix,X,as.logical(zeroDiag),as.integer(threads))

#' Directional normalised mutual information matrix
#'
#' Calculates directed normalised mutual information between each two features, that is
#' \deqn{\frac{I(X_i,X_j)}{H(X_j)}.}
#' @template input
#' @template matrix
#' @examples
#' dnmiMatrix(iris)
#' @export
dnmiMatrix<-function(X,zeroDiag=TRUE,threads=0)
 .Call(C_dnmiMatrix,X,as.logical(zeroDiag),as.integer(threads))

#' Conditional mutual information scores
#'
#' Calculates conditional mutual information between each features and the decision, that is
#' \deqn{I(X;Y|Z).}
#' @template input
#' @template y
#' @param Z Condition; should be given as a factor, but other options are accepted, as for features. 
#' @return A numerical vector with conditional mutual information scores, with names copied from \code{X}.
#' @examples
#' cmiScores(iris[,-5],iris$Species,iris$Sepal.Length)
#' @export
cmiScores<-function(X,Y,Z,threads=0)
 .Call(C_cmi,X,Y,Z,as.integer(threads))

#' Conditional mutual information matrix
#'
#' Calculates conditional mutual information between each two features given another one, that is
#' \deqn{I(X_i;X_j|Z).}
#' @template input
#' @param Z Condition; should be given as a factor, but other options are accepted, as for features. 
#' @template matrix
#' @examples
#' cmiMatrix(iris[,-5],iris[,5])
#' @export
cmiMatrix<-function(X,Z,zeroDiag=TRUE,threads=0)
 .Call(C_cmiMatrix,X,Z,as.logical(zeroDiag),as.integer(threads))

#' Joint mutual information scores
#'
#' Calculates joint mutual information between each feature joint with some other vector \code{Z} with the decision, that is
#' \deqn{I(X,Z;Y).}
#' This is the same as conditional mutual information between X and Y plus a constant that depends on Y and Z, that is
#' \deqn{I(X,Z;Y)=I(X;Y|Z)+I(Y;Z).}
#' @template input
#' @template y
#' @param Z Other vector; should be given as a factor, but other options are accepted, as for features. 
#' @return A numerical vector with joint mutual information scores, with names copied from \code{X}.
#' @examples
#' jmiScores(iris[,-5],iris$Species,iris$Sepal.Length)
#' @export
jmiScores<-function(X,Y,Z,threads=0)
 .Call(C_jmi,X,Y,Z,as.integer(threads))

#' Joint mutual information matrix
#'
#' Calculates mutual information between each feature and a joint mix of each other feature with a given feature, that is
#' \deqn{I(X_i;X_j,Z).}
#' @template input
#' @param Z Condition; should be given as a factor, but other options are accepted, as for features. 
#' @template matrix
#' @examples
#' jmiMatrix(iris[,-5],iris[,5])
#' @export
jmiMatrix<-function(X,Z,zeroDiag=TRUE,threads=0)
 .Call(C_jmiMatrix,X,Z,as.logical(zeroDiag),as.integer(threads))

#' Normalised joint mutual information scores
#'
#' Calculated normalised mutual information between each feature joint with some other vector \code{Z} and the decision, that is
#' \deqn{\frac{I(X,Z;Y)}{H(X,Y,Z)}.}
#' This is the same as in the criterion used by \code{\link{DISR}} and \code{\link{NJMIM}}.
#' @template input
#' @template y
#' @param Z Other vector; should be given as a factor, but other options are accepted, as for features. 
#' @return A numerical vector with the normalised joint mutual information scores, with names copied from \code{X}.
#' @examples
#' njmiScores(iris[,-5],iris$Species,iris$Sepal.Length)
#' @export
njmiScores<-function(X,Y,Z,threads=0)
 .Call(C_njmi,X,Y,Z,as.integer(threads))

#' Normalised joint mutual information matrix
#'
#' Calculates normalised mutual information between each feature and a joint mix of each other feature with a given feature, that is
#' \deqn{\frac{I(X_i;X_j,Z)}{H(X_i,X_j,Z)}.}
#' @template input
#' @param Z Condition; should be given as a factor, but other options are accepted, as for features. 
#' @template matrix
#' @examples
#' njmiMatrix(iris[,-5],iris[,5])
#' @export
njmiMatrix<-function(X,Z,zeroDiag=TRUE,threads=0)
 .Call(C_njmiMatrix,X,Z,as.logical(zeroDiag),as.integer(threads))

#' Gini impurity scores
#'
#' Calculates Gini impurity between each feature and the decision, that is
#' \deqn{G(X;Y)=\sum_{xy} \frac{p_{xy}^2}{p_x}-\sum_y p_y^2.}
#' @template input
#' @template y
#' @return A numerical vector with Gini impurity scores, with names copied from \code{X}.
#' @examples
#' impScores(iris[,-5],iris$Species)
#' @export
impScores<-function(X,Y,threads=0)
 .Call(C_im,X,Y,as.integer(threads))

#' Entropy scores
#'
#' Calculates entropy of each feature, that is
#' \deqn{H(X).}
#' @template input
#' @return A numerical vector with entropy scores, with names copied from \code{X}.
#' @examples
#' hScores(iris[,-5])
#' @export
hScores<-function(X,threads=0)
 .Call(C_h,X,as.integer(threads))

#' Joint entropy scores
#'
#' Calculates joint entropy of each feature and a condition \code{Y}, that is
#' \deqn{H(X,Y).}
#' @template input
#' @template y
#' @return A numerical vector with entropy scores, with names copied from \code{X}.
#' @examples
#' jhScores(iris[,-5],iris[,5])
#' @export
jhScores<-function(X,Y,threads=0)
 .Call(C_jh,X,Y,as.integer(threads))

#' Mutual information of feature triples
#'
#' Calculates mutual information of each triple of features, that is
#' \deqn{I(X_i;X_j;X_k).}
#' @template input
#' @return A data frame with four columns; first three (\code{Var1}, \code{Var2} and \code{Var3}) are names of features, fourth, \code{MI} is the value of the mutual information.
#' The order of features does not matter, hence only \deqn{n(n-1)(n-2)/6} unique, sorted triples are evaluated.
#' @note In a current version, the maximal number of features accepted is 2345, which gives a bit less than 2^32 triples.
#' The equation used for calculation is
#' \deqn{I(X_i;X_j;X_k)=I(X_i;X_k)+I(X_j;X_k)-I(X_i,X_j;X_k).}
#' Henceforth, please mind that rounding errors may occur and influence reproducibility.
#' @examples
#' triScores(iris)
#' @export
triScores<-function(X,threads=0)
 data.frame(.Call(C_tri,X,as.integer(threads)))
