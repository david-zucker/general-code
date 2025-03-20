# EXAMPLE OF EDF AND LILLIEFORS'S TEST

#FUNCTION TO COMPUTE THE EMPIRICAL DISTRIBUTION FUNCTION OF A DATA SAMPLE
my.edf = function(x) {
  n = length(x)
  z = sort(unique(x))
  m = length(z)
  f1 = NULL
  f2 = NULL
  for (j in 1:m) {
    z.cur = z[j]
    f1.cur = sum(x < z.cur)/n
    f2.cur = sum(x <= z.cur)/n
    f1 = c(f1,f1.cur)
    f2 = c(f2,f2.cur)
  }
  res = list(z=z, f1=f1, f2=f2)
}

#FUNCTION TO COMPUTE LILLIEFORS STATISTIC 
my.lilli = function(x) {
  meanval = mean(x)
  sdval = sd(x)
  edf.res = my.edf(x)
  prb = pnorm(edf.res$z,meanval,sdval)
  d1 = max(abs(edf.res$f1-prb))
  d2 = max(abs(edf.res$f2-prb))
  d = max(c(d1,d2))
  return(d)
}

#FUNCTION TO COMPUTE THE P-VALUE FOR A GIVEN VALUE OF THE LILLIEFORS STATISTIC
#VIA SIMULATION 
lilli.pval = function(d.obs,n,nsim) {
  llsim = NULL
  for (s in 1:nsim) {
    x = rnorm(n,0,1)
    d = my.lilli(x)
    llsim = c(llsim,d)
  }
  pval = mean(llsim >= d.obs)
  return(pval)
}  

#EXAMPLE

library(nortest)

xsam = rgamma(20,2,1)
  
#USING MY CODE
dd = my.lilli(xsam)
pval = lilli.pval(dd,20,100000)
print(cbind(dd,pval))

#USING R FUNCTION FROM nortest PACKAGE
print(lillie.test(xsam))
