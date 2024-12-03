#PROGRAM TO DEMONSTRATE NONPARAMETRIC BOOTSTRAP
#PARAMETER OF INTEREST IS COEFFICIENT OF VARIATION SD(Y)/E[Y]

#PRELIMINARIES
set.seed(342345)
n = 50
alftru = 3
bettru = 0.5
M = 1000

#TRUE PARAMETER VALUE is 1/sqrt(alftru) = 0.5774

#GENERATE SAMPLE OF SIZE n FROM GAMMA(alftru,bettru) DISTRIBUTION
#using parameterization in Wasserman
yorig = rgamma(n,alftru,scale=bettru)

#COMPUTE ESTIMATES OF THETA FROM THE ORIGINAL DATA
thhat = sd(yorig)/mean(yorig)

#BOOTSTRAP
tstar = rep(0,M)
for (m in 1:M) {
  ycur = sample(yorig, size=n, replace=TRUE) 
  thhst = sd(ycur)/mean(ycur)
  tstar[m] = thhst - thhat
 }
q = quantile(tstar,c(0.975,0.025),type=4)
bootci = thhat - q
print(cbind(bootci))

###########################################################

#SAMPLE OUTPUT

#bootci
#97.5% 0.4219416
#2.5%  0.6215228


