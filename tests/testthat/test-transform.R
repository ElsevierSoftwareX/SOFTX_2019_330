context("test-transform.R")

test_that("Kendall transformation behaves properly",{
 for(n in 5:10){
  x<-runif(n)
  tx<-kTransform(x)
  expect_equal(tx,kTransform(x*2))
  expect_equal(length(tx),n*(n-1))

  cx<-cut(x,2)
  kTransform(cx)->tcx
  expect_true(is.factor(tcx))
  expect_equal(length(tcx),n*(n-1))
  expect_error(
   kTransform(cut(x,3)),
   "Unordered factor with more"
  )
  expect_equal(levels(tcx),c("<",">","="))

  kTransform(1:n)->kn
  kTransform(n:1)->kr
  expect_identical(sort(kn),sort(kr))
  expect_equal(length(table(factor(kn))),2)
  kTransform(rep(pi,n))->kz
  expect_true(all(kz==kz[1]))
 }
})

test_that("Kendall transformation handles various inputs",{
  expect_error(
   kTransform("abc"),
   "Unsupported"
  )
  expect_error(
   kTransform(data.frame(a=letters,stringsAsFactors=FALSE)),
   "Unsupported"
  )
  expect_error(
   kTransform(NULL),
   "Unsupported"
  )
  eo4<-factor(c("<","<","<",">","=","<",">","=","<",">",">",">"),levels=c("<",">","="))
  eo3<-factor(c("<","<",">","<",">",">"),levels=c("<",">","="))
  eo2<-factor(c("<",">"),levels=c("<",">","="))
  eona<-factor(c(NA,"<",NA,NA,">",NA),levels=c("<",">","="))
  expect_equal(kTransform(c(FALSE,TRUE)),eo2)
  expect_equal(kTransform(c(-Inf,0,Inf)),eo3)
  expect_equal(kTransform(c(1.1,2.2,2.2,3.3)),eo4)
  expect_equal(kTransform(c(1.1,NA,3.3)),eona)
  expect_equal(kTransform(c(1.1,NaN,3.3)),eona)
  expect_equal(kTransform(c(FALSE,NA,TRUE)),eona)
  expect_equal(kTransform(c(0,NA,1)),eona)
  expect_equal(kTransform(ordered(c("a","b","b","c"))),eo4)
})

test_that("Kendall transformation sanity test",{
 expect_equal(
  sort(as.character(kTransform(c(1,1,2,0)))),
  sort(c('=','>','<','=','>','<','<','<','<','>','>','>'))
 )
})

test_that("Kendall works on data frames",{
 x<-iris[,-5]
 expect_equal(
  kTransform(x),
  data.frame(lapply(x,kTransform))
 )
})

test_that("Kendall transformation is reversible",{
 for(e in 1:20){
  set.seed(e)
  x<-sample(runif(1,1,20),runif(1,3,20),replace=TRUE)
  kx<-kTransform(x)
  expect_equal(kInverse(kx),rank(x))
 }
})

test_that("Kendall transformation inverses score vectors",{
 expect_equal(
  kInverse((kTransform(1:5)==">")*pi),
  1:5
 )
})

test_that("Inverse KT throws proper errors",{
 expect_error(kInverse(iris$Species),"Factor does not seem")
 expect_error(kInverse(1:5),"Invalid size")
 expect_error(kInverse(iris),"Invalid input")
 expect_error(kInverse(c(1i,2i)),"Invalid input")
 expect_error(kInverse(integer(0)),"Input too short")
 expect_error(kInverse(1),"Input too short")
})

