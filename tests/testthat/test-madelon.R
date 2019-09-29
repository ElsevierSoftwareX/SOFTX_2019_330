context("test-madelon.R")

data(MadelonD)

set.seed(17)

input<-list(X=MadelonD$X,Y=MadelonD$Y,k=20)
pmr<-sample(nrow(MadelonD$X),nrow(MadelonD$X))
pmc<-sample(ncol(MadelonD$X),ncol(MadelonD$X))
inputp<-list(X=MadelonD$X[pmr,pmc],Y=MadelonD$Y[pmr],k=20)

for(algo in c("MIM","JMIM","NJMIM","JMI","DISR","CMIM","MRMR","JIM")){
 test_that(sprintf("Native %s works the same on permuted Madelon",algo,algo),{
  do.call(algo,input)->j
  do.call(algo,inputp)->p
  expect_equal(j$score,p$score)
  expect_equal(names(j$selection),names(p$selection))
 })
}

# We need exception for CMI since attributes become indistinguishable for it
test_that(sprintf("Native %s works the same on permuted Madelon",algo,algo),{
 inputp<-list(X=MadelonD$X[pmr,],Y=MadelonD$Y[pmr],k=20)
 do.call(CMI,input)->j
 do.call(CMI,inputp)->p
 expect_equal(j$score,p$score)
 expect_equal(names(j$selection),names(p$selection))
})
