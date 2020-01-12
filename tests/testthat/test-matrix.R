context("test-matrix.R")

pureMiMatrix<-function(X,zeroDiag=TRUE){
 sapply(1:ncol(X),function(e) miScores(X,X[,e]))->ans
 if(zeroDiag) diag(ans)<-0
 colnames(ans)<-names(X)
 ans
}

pureNmiMatrix<-function(X,zeroDiag=TRUE){
 sapply(1:ncol(X),function(e) njmiScores(X,X[,e],rep(1,nrow(X))))->ans
 if(zeroDiag) diag(ans)<-0
 colnames(ans)<-names(X)
 ans
}

pureDnmiMatrix<-function(X,zeroDiag=TRUE){
 sapply(1:ncol(X),function(e) miScores(X,X[,e]))->ans
 if(zeroDiag) diag(ans)<-0
 colnames(ans)<-names(X)
 t(t(ans)/hScores(X))
}

pureCmiMatrix<-function(X,Z,zeroDiag=TRUE){
 sapply(1:ncol(X),function(e) cmiScores(X,X[,e],Z))->ans
 if(zeroDiag) diag(ans)<-0
 colnames(ans)<-names(X)
 ans
}

pureJmiMatrix<-function(X,Z,zeroDiag=TRUE){
 sapply(1:ncol(X),function(e) jmiScores(X,X[,e],Z))->ans
 if(zeroDiag) diag(ans)<-0
 colnames(ans)<-names(X)
 ans
}

pureNjmiMatrix<-function(X,Z,zeroDiag=TRUE){
 sapply(1:ncol(X),function(e) njmiScores(X,X[,e],Z))->ans
 if(zeroDiag) diag(ans)<-0
 colnames(ans)<-names(X)
 ans
}

test_that("mi matrix works",{
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

test_that("nmi matrix works as pure",{
 expect_equal(
  pureNmiMatrix(iris,FALSE),
  nmiMatrix(iris,FALSE)
 )
})

test_that("dnmi matrix works as pure",{
 expect_equal(
  pureDnmiMatrix(iris,FALSE),
  dnmiMatrix(iris,FALSE)
 )
})

test_that("cmi matrix works as pure",{
 expect_equal(
  pureCmiMatrix(iris,iris[,5],FALSE),
  cmiMatrix(iris,iris[,5],FALSE)
 )
})

test_that("jmi matrix works as pure",{
 expect_equal(
  pureJmiMatrix(iris,iris[,5],FALSE),
  jmiMatrix(iris,iris[,5],FALSE)
 )
})

test_that("njmi matrix works as pure",{
 expect_equal(
  pureNjmiMatrix(iris,iris[,5],FALSE),
  njmiMatrix(iris,iris[,5],FALSE)
 )
})

test_that("triScores works",{
 triScores(cbind(iris,C=rep(1,150)))->z
 expect_equal(names(z),c("Var1","Var2","Var3","MI"))
 expect_equal(nrow(z),6*5*4/6)

 z[z$Var3=="C",]->z
 expect_equal(z$MI,rep(0,10))

 expect_equal(triScores(iris[,c(1,2,1)])$MI,miScores(iris[,1],iris[,2]))
})
