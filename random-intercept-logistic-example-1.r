#TOENAIL DATA
library(GLMMadaptive)
dat1 = read.table('toenaildat.txt',header=T)
print(dat1[1:10,])
dat1$trt.time = dat1$treatn * dat1$time
res1 = mixed_model(y ~ time + treatn + trt.time, ~ 1|idnew, data=dat1, family=binomial(),
  control=list(nAGQ=20))
print(summary(res1))
Omega1 = solve(res1$Hessian)
print(Omega1)