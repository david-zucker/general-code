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
 
#GRADIENT OF LOG-LIKELIHOOD FUNCTION
llgrad = function(phi) {
  alf = phi[1]
  bet = phi[2]
  r1 = mean(log(ycur))
  r2 = mean(ycur)
  d1a = -log(bet) - digamma(alf) + r1
  d1b = -(alf/bet) + r2/(bet^2)
  llgrad = c(d1a,d1b)
} 

#ASYMPTOTIC VARIANCE CALCULATION
asyvar = function(alfhat,bethat) {
  infmat = matrix(0,2,2)
  infmat[1,1] = trigamma(alfhat)
  infmat[1,2] = 1/bethat
  infmat[2,1] = infmat[1,2]
  infmat[2,2] = alfhat/(bethat^2)
  vmat = solve(infmat)
  grd = matrix(c(bethat,alfhat),2,1)
  asyvar = t(grd) %*% vmat %*% grd / n
  if (asyvar < 0) browser()
  return(asyvar)
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
  method="L-BFGS-B", control=list(fnscale=-1),hessian=T)
alfhat = mle$par[1]
bethat = mle$par[2]
thhat = alfhat*bethat

#ASYMPTOTIC CI WITH ASYMPTOTIC VARIANCE
asyvar_orig = asyvar(alfhat,bethat)
asysd = sqrt(asyvar_orig)
cihw1 = za*asysd
theta_lo_1 = thhat - cihw1
theta_hi_1 = thhat + cihw1

#BOOTSTRAP
tstar = rep(0,B)
wstar = rep(0,B)
for (b in 1:B) {
  ycur = rgamma(n,alfhat,scale=bethat)
  ybar = mean(ycur)
  yvar = var(ycur)
  bet0 = yvar/ybar
  alf0 = ybar/bet0
  mle = optim(c(alf0,bet0),llik,gr=llgrad,lower=c(1e-4,1e-4),
    method="L-BFGS-B", control=list(fnscale=-1), hessian=T)
  alfhst = mle$par[1]
  bethst = mle$par[2]
  thhst = alfhst*bethst
  asyvar_cur = asyvar(alfhst,bethst)
  tstar[b] = thhst - thhat
  wstar[b] = (thhst - thhat) / sqrt(asyvar_cur)
}

#ASYMPTOTIC CI WITH BOOTSTRAP VARIANCE
cihw2 = za*sd(tstar)
theta_lo_2 = thhat - cihw2
theta_hi_2 = thhat + cihw2

#BOOSTRAP CI
q = quantile(tstar,c(0.975,0.025),type=4)
theta_lo_3 = thhat - q[1]
theta_hi_3 = thhat - q[2]

#STUDENTIZED BOOSTRAP CI
q = quantile(wstar,c(0.975,0.025),type=4)
theta_lo_4 = thhat - q[1]*asysd
theta_hi_4 = thhat - q[2]*asysd

#SUMMARIZE RESULTS
rslt = rbind(
  c(theta_lo_1,theta_hi_1),
  c(theta_lo_2,theta_hi_2),
  c(theta_lo_3,theta_hi_3),
  c(theta_lo_4,theta_hi_4))
rslt = as.data.frame(rslt)
colnames(rslt) = c('theta_lo','theta_hi')
rownames(rslt) = c('asy/asy', 'asy/boot', 'bootstrap', 'stu boot')

cat('\n')
print(cbind(thtru,thhat))
cat('\n')
print(rslt)
cat('\n')



###########################################################

#SAMPLE OUTPUT


#Sample size of n=50

#  thtru    thhat
#      2 2.270189

#          theta_lo theta_hi
#asy/asy   1.866296 2.674081
#asy/boot  1.876934 2.663444
#bootstrap 1.851044 2.643153
#stu boot  1.906108 2.744996


#Sample size of n=200

#  thtru    thhat
#      2 1.940273

#          theta_lo theta_hi
#asy/asy   1.748640 2.131906
#asy/boot  1.758756 2.121791
#bootstrap 1.750191 2.119082
#stu boot  1.761410 2.143114


