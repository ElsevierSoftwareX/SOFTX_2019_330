context("test-tools.R")

test_that("joinf",{
 joinf(
  c(T,F,T,F,F),
  data.frame(
   i=factor(c("a","a","b","b","b")),
   j=rep(17,5)
  )
 )->z
 expect_equal(z,factor(sprintf("l%d",c(1:4,4))))
})

