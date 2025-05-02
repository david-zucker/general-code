#PROGRAM TO DEMONSTRATE STUDENTIZED BOOTSTRAP 
#PARAMETER OF INTEREST IS E[Y]

#gamma distribution with Wasserman parameterization

#PRELIMINARIES
set.seed(342345)
za = -qnorm(0.025)
n = 30
alftru = 3
bettru = 0.5
M = 200

#TRUE PARAMETER VALUE IS 1.5

#GENERATE SAMPLE OF SIZE n FROM GAMMA(alftru,bettru) DISTRIBUTION
yorig = rgamma(n,alftru,scale=bettru)

#COMPUTE ESTIMATES OF THETA FROM THE ORIGINAL DATA
thhat = mean(yorig)

#VARIANCE ESTIMATE
estsd = sd(yorig)/sqrt(n)

#BOOTSTRAP LOOP
tstar = rep(0,M)
for (m in 1:M) {
  ycur = sample(yorig, size=n, replace=TRUE) 
  thhst = mean(ycur)
  estsd_boot = sd(ycur)/sqrt(n)
  tstar[m] = (thhst - thhat)/estsd_boot
}

#COMPUTE CI
q = quantile(tstar,c(0.975,0.025),type=4)
thlo = thhat - q[1]*estsd
thhi = thhat - q[2]*estsd
ans = cbind(thhat,thlo,thhi)
print(ans)

###########################################################

#SAMPLE OUTPUT

#     thhat      thlo      thhi
#  1.584345  1.345511  1.902669

