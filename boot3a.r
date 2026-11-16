#PROGRAM TO DEMONSTRATE PAIRED-SAMPLE NONPARAMETRIC BOOTSTRAP

#PRELIMINARIES
library(MASS)
set.seed(342345)
za = -qnorm(0.025)
n = 45
rhotru = 0.7
M = 1000

#GENERATE ORIGINAL DATA AND COMPUTE ESTIMATE
mu = c(0,0)
Sigma = matrix(c(1,rhotru,rhotru,1),2,2)
yorig = mvrnorm(n,mu,Sigma) 
y1orig = yorig[,1]
y2orig = yorig[,2]
rhohat = cor(y1orig,y2orig)

#BOOTSTRAP
tstar = rep(0,M)
for (m in 1:M) {
  idxcur = sample.int(n,size=n,replace=TRUE)
  y1 = yorig[idxcur,1]
  y2 = yorig[idxcur,2]
  rhhst = cor(y1,y2)
  tstar[m] = rhhst - rhohat
 }
q = quantile(tstar,c(0.975,0.025),type=4)
bootci = rhohat - q
bootci = matrix(bootci,1,2)
bootci = as.data.frame(bootci)
colnames(bootci) = c('rho_lo','rho_hi')

#PRINT OUTPUT
print(cbind(rhotru,rhohat))
cat('\n')
print(bootci)
###########################################################

#SAMPLE OUTPUT
# 
#       rhotru   rhohat
#  [1,]    0.7 0.624156
#
#       rho_lo    rho_hi
#  1 0.4330863 0.9290806