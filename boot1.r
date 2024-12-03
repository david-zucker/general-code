#PROGRAM TO DEMONSTRATE ASYMPTOTIC CI WITH BOOTSTRAP VARIANCE
#PARAMETER OF INTEREST IS COEFFICIENT OF VARIATION SD(Y)/E[Y]

#PRELIMINARIES
set.seed(342345)
za = -qnorm(0.025)
n = 50
alftru = 3
bettru = 2
M = 200

#TRUE PARAMETER VALUE is 1/sqrt(alftru) = 0.5774

#GENERATE SAMPLE OF SIZE n FROM GAMMA(alftru,bettru) DISTRIBUTION
yorig = rgamma(n,alftru,rate=bettru)

#COMPUTE ESTIMATES OF THETA FROM THE ORIGINAL DATA
thhat = sd(yorig)/mean(yorig)

#BOOTSTRAP
tstar = rep(0,M)
for (m in 1:M) {
  ycur = sample(yorig, size=n, replace=TRUE) 
  thhst = sd(ycur)/mean(ycur)
  tstar[m] = thhst
 }
sdhat = sd(tstar)

#ASYMPTOTIC CI WITH BOOTSTRAP VARIANCE
cihw = za*sdhat
thlo = thhat - cihw
thhi = thhat + cihw
ans = cbind(thhat,thlo,thhi)
print(ans)

###########################################################

#SAMPLE OUTPUT

#     thhat      thlo      thhi
# 0.5212792 0.4273496 0.6152089

