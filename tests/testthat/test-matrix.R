context("test-matrix.R")

test_that("pure mi matrix works",{
 pureMiMatrix(iris,FALSE)->M
 expect_equal(diag(M),hScores(iris))
 expect_equal(t(M),M)
 expect_equal(rownames(M),colnames(M))
 expect_equal(rownames(M),names(iris))
 expect_equal(
  diag(pureMiMatrix(iris)),
  setNames(rep(0,5),names(iris))
 )
})
