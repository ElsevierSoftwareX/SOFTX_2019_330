context("test-ensemble.R")

test_that("EJMI stability",{
 set.seed(7)
 EJMI(iris[,-5],iris[,5],k=3,n=20,p=.5,th=1)->a
 set.seed(7)
 EJMI(iris[,-5],iris[,5],k=3,n=20,p=.5,th=0)->b
 expect_equal(a,b)

 set.seed(77)
 EJMI(iris[,-5],iris[,5],k=3,n=20,p=.5,th=0)->bp
 expect_false(identical(bp,b))
})
