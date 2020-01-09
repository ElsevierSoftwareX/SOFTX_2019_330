context("test-transform.R")

test_that("Kendall transform behaves properly",{
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

test_that("Kendall transform handles various inputs",{
  expect_error(
   kTransform("abc"),
   "Unsupported"
  )
  expect_error(
   kTransform(NULL),
   "Unsupported"
  )
  eo3<-factor(c("<","<","<",">","=","<",">","=","<",">",">",">"),levels=c("<",">","="))
  eo2<-factor(c("<",">"),levels=c("<",">","="))
  expect_equal(kTransform(c(FALSE,TRUE)),eo2)
  expect_equal(kTransform(c(1.1,2.2,2.2,3.3)),eo3)
  expect_equal(kTransform(ordered(c("a","b","b","c"))),eo3)
})

test_that("Kendall transform sanity test",{
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

