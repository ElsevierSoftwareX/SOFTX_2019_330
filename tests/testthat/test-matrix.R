context("test-matrix.R")

pureMiMatrix<-function(X,zeroDiag=TRUE){
 sapply(1:ncol(X),function(e) miScores(X,X[,e]))->ans
 if(zeroDiag) diag(ans)<-0
 colnames(ans)<-names(X)
 ans
}

test_that("Mi matrix works",{
 miMatrix(iris,FALSE)->M
 expect_equal(diag(M),hScores(iris))
 expect_equal(t(M),M)
 expect_equal(rownames(M),colnames(M))
 expect_equal(rownames(M),names(iris))
 expect_equal(
  diag(miMatrix(iris)),
  setNames(rep(0,5),names(iris))
 )
})

test_that("mi matrix works as pure",{
 expect_equal(
  pureMiMatrix(iris,FALSE),
  miMatrix(iris,FALSE)
 )
})
