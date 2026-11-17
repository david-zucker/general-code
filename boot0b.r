#PROGRAM TO DEMONSTRATE SIMULATED STANDARD DEVIATION
#PARAMETER OF INTEREST IS LOG VARIANCE

#gamma distribution with Wasserman parameterization

#PRELIMINARIES
set.seed(342345)
n = 100
alf = 3
bet = 0.5
R = 1000

#ORIGINAL DATA
yorig = rgamma(n,alf,scale=bet)

#We act as if we don't know the distribution from which the original data came,
#which is the usual situation in practice.

#BOOTSTRAP
lthhat = rep(0,R)
for (r in 1:R) {
  ycur = sample(yorig,size=n,replace=TRUE)
  lthh = log(var(ycur))
  lthhat[r] = lthh
 }
sd_log_th_sim = sd(lthhat)

print(sd_log_th_sim)

