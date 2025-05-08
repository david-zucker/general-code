#PROGRAM TO DEMONSTRATE PARAMETRIC BOOTSTRAP
#GAMMA-DISTRIBUTED DATA
#PARAMETER OF INTEREST IS THE EXPECTATION alpha*beta

#PRELIMINARIES
set.seed(342345)
za = -qnorm(0.025)
n = 50
alftru = 2
bettru = 1
thtru = alftru*bettru
B = 1000

#LOG-LIKELIHOOD FUNCTION
llik = function(phi) {
  alf = phi[1]
  bet = phi[2]
  r1 = mean(log(ycur))
  r2 = mean(ycur)
  llik = -alf * log(bet) - lgamma(alf) + (alf-1)*r1 - r2/bet
 }
 
llgrad = function(phi) {
  alf = phi[1]
  bet = phi[2]
  r1 = mean(log(ycur))
  r2 = mean(ycur)
  d1a = log(bet) - digamma(alf) + r1
  d1b = -(alf/bet) + r2/(bet^2)
  llgrad = c(d1a,d1b)
 } 
  
#GENERATE SAMPLE OF SIZE n FROM GAMMA(alftru,bettru) DISTRIBUTION
yorig = rgamma(n,alftru,scale=bettru)

#COMPUTE ESTIMATES OF ALPHA, BETA, AND THETA FROM THE ORIGINAL DATA
ycur = yorig
ybar = mean(ycur)
yvar = var(ycur)
bet0 = yvar/ybar
alf0 = ybar/bet0
mle = optim(c(alf0,bet0),llik,gr=llgrad,lower=c(1e-4,1e-4),
  method="L-BFGS-B", control=list(fnscale=-1))
alfhat = mle$par[1]
bethat = mle$par[2]
thhat = alfhat*bethat

#ASYMPTOTIC CI WITH ASYMPTOTIC VARIANCE
infmat = matrix(0,2,2)
infmat[1,1] = trigamma(alfhat)
infmat[1,2] = 1/bethat
infmat[2,1] = infmat[1,2]
infmat[2,2] = alfhat/(bethat*2)
vmat = solve(infmat)
grd = matrix(c(bethat,alfhat),2,1)
asyvar = t(grd) %*% vmat %*% grd / n
cihw1 = za*sqrt(asyvar)
theta_lo_1 = thhat - cihw1
theta_hi_1 = thhat + cihw1

#BOOTSTRAP
tstar = rep(0,B)
for (b in 1:B) {
  ycur = rgamma(n,alfhat,scale=bethat)
  ybar = mean(ycur)
  yvar = var(ycur)
  bet0 = yvar/ybar
  alf0 = ybar/bet0
  mle = optim(c(alf0,bet0),llik,gr=llgrad,lower=c(1e-4,1e-4),
    method="L-BFGS-B", control=list(fnscale=-1))
  alfhst = mle$par[1]
  bethst = mle$par[2]
  thhst = alfhst*bethst
  tstar[b] = thhst - thhat
 }

#ASYMPTOTIC CI WITH BOOTSTRAP VARIANCE
cihw2 = za*sd(tstar)
theta_lo_2 = thhat - cihw2
theta_hi_2 = thhat + cihw2

#BOOSTRAP CI
q = quantile(tstar,c(0.975,0.025),type=4)
theta_lo_3 = thhat - q[1]
theta_hi_3 = thhat - q[2]

#SUMMARIZE RESULTS
rslt = rbind(
  c(theta_lo_1,theta_hi_1),
  c(theta_lo_2,theta_hi_2),
  c(theta_lo_3,theta_hi_3))
rslt = as.data.frame(rslt)
colnames(rslt) = c('theta_lo','theta_hi')
rownames(rslt) = c('asy/asy', 'asy/boot', 'bootstrap')

cat('\n')
print(cbind(thtru,thhat))
cat('\n')
print(rslt)
cat('\n')



###########################################################

#SAMPLE OUTPUT


#Sample size of n=50

#     thtru    thhat
#[1,]     2 2.278366
#
#          theta_lo theta_hi
#asy/asy   1.936886 2.619846
#asy/boot  1.820196 2.736537
#bootstrap 1.699417 2.587763


#Sample size of n=200

#     thtru    thhat
#[1,]     2 1.943921
#
#          theta_lo theta_hi
#asy/asy   1.797398 2.090445
#asy/boot  1.739997 2.147846
#bootstrap 1.717387 2.119066




