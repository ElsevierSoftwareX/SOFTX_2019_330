devtools::load_all()
library(microbenchmark)
data(MadelonD)
with(MadelonD,cbind(cmiScores2(X,X[,23],Y),cmiScores(X,X[,23],Y)))->ans
stopifnot(all.equal(ans[,1],ans[,2]))
print(with(MadelonD,microbenchmark(cmiScores2(X,X[,23],Y),cmiScores(X,X[,23],Y,1),times=10)))
print(microbenchmark(cmiScores2(iris,iris[,2],iris[,5]),cmiScores(iris,iris[,2],iris[,5]),times=100))
