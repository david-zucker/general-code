#PROGRAM TO DEMONSTRATE ASYMPTOTIC CI WITH BOOTSTRAP VARIANCE

#PRELIMINARIES
library(MASS)
set.seed(342345)
za = -qnorm(0.025)
n = 20
nreps = 1000
rhotru = 0.7
M = 100

#PREPARE FOR SIMULATION
cover = 0
lngth = NULL

#SIMULATION LOOP
for (irep in 1:nreps) {

#GENERATE ORIGINAL DATA AND COMPUTE ESTIMATE
mu = c(0,0)
Sigma = matrix(c(1,rhotru,rhotru,1),2,2)
yorig = mvrnorm(n,mu,Sigma)/rchisq(1,2) 
y1orig = yorig[,1]
y2orig = yorig[,2]
rhohat = cor(y1orig,y2orig)

#BOOTSTRAP VARIANCE
tstar = rep(0,M)
for (m in 1:M) {
  idxcur = sample.int(n,size=n,replace=TRUE)
  y1 = yorig[idxcur,1]
  y2 = yorig[idxcur,2]
  rhhst = cor(y1,y2)
  tstar[m] = rhhst - rhohat
 }
sdhat = sd(tstar)

#ASYMPTOTIC CI WITH BOOTSTRAP VARIANCE
cihw = za*sdhat
thlo = rhohat - cihw
thhi = rhohat + cihw
ci = c(thlo,thhi)
cover = cover + (ci[1] <= fshtru)*(fshtru <= ci[2])
cilng = 2*cihw
lngth = c(lngth,cilng)

}
#END OF SIMULATION LOOP

#SUMMARY RESULTS
cover = cover/nreps
lmean = mean(lngth)
lsd = sd(lngth)
print(cbind(n,M,nreps,cover,lmean,lsd))

###########################################################

#SAMPLE OUTPUT 
#      n   M nreps cover     lmean       lsd
#[1,] 20 100  1000 0.926 0.9217707 0.1764025
