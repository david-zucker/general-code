#PROGRAM TO DEMONSTRATE PAIRED-SAMPLE NONPARAMETRIC BOOTSTRAP

#PRELIMINARIES
library(MASS)
set.seed(342345)
za = -qnorm(0.025)
n = 20
rhotru = 0.7
M = 1000

fishz = function(r) {
  fishz = 0.5 * (log(1+r) - log(1-r))
}

fishz.inv = function(f) {
  q = exp(2*f)
  ans = (q-1)/(q+1)
  return(ans)
}

#GENERATE ORIGINAL DATA AND COMPUTE ESTIMATE
fshtru = fishz(rhotru)
mu = c(0,0)
Sigma = matrix(c(1,rhotru,rhotru,1),2,2)
yorig = mvrnorm(n,mu,Sigma) 
y1orig = yorig[,1]
y2orig = yorig[,2]
rhohat = cor(y1orig,y2orig)
fshhat = fishz(rhohat)

#ASYMPTOTIC CI
cihw = za/sqrt(n)
thlo = fshhat - cihw
thhi = fshhat + cihw
asyci.fishz = c(thlo,thhi)
asyci.rho = fishz.inv(asyci.fishz)

#BOOTSTRAP
tstar = rep(0,M)
for (m in 1:M) {
  idxcur = sample.int(n,size=n,replace=TRUE)
  y1 = yorig[idxcur,1]
  y2 = yorig[idxcur,2]
  rhhst = cor(y1,y2)
  fshhst = fishz(rhhst)
  tstar[m] = fshhst - fshhat
 }
q = quantile(tstar,c(0.975,0.025),type=4)
bootci.fishz = fshhat - q 
bootci.rho = fishz.inv(bootci.fishz)

#PRINT OUTPUT
print(cbind(rhotru,rhohat))
cat('\n')
print(cbind(fshtru,fshhat))
cat('\n')
print(cbind(asyci.rho,bootci.rho))
cat('\n')
print(cbind(asyci.fishz,bootci.fishz))
cat('\n')

###########################################################

#SAMPLE OUTPUT
# 
#      rhotru    rhohat
# [1,]    0.7 0.5386438
# 
#         fshtru    fshhat
# [1,] 0.8673005 0.6022431
# 
#       asyci.rho   bootci.rho
# 97.5% 0.1625276  -0.08377435
#  2.5%  0.7780872  0.79649247
# 
#       asyci.fishz bootci.fishz
# 97.5%   0.1639818  -0.08397116
# 2.5%    1.0405044   1.08894419






