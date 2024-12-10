#PROGRAM TO DEMONSTRATE TWO-SAMPLE PARAMETRIC STUDENTIZED BOOTSTRAP

#PRELIMINARIES
set.seed(342345)
za = -qnorm(0.025)
n = 20
alftru1 = 2
bettru1 = 1
exptd1 = alftru1/bettru1
alftru2 = 0.5
bettru2 = 1
exptd2 = alftru2/bettru2
thtru = exptd1/exptd2
M = 1000

#LOG-LIKELIHOOD FUNCTION
llik = function(phi) {
  alf = phi[1]
  bet = phi[2]
  r1 = mean(log(ycur))
  r2 = mean(ycur)
  ll = alf * log(bet) - lgamma(alf) + (alf-1)*r1 - bet*r2
  llik = -ll
 }
 
llgrad = function(phi) {
  alf = phi[1]
  bet = phi[2]
  r1 = mean(log(ycur))
  r2 = mean(ycur)
  d1a = log(bet) - digamma(alf) + r1
  d1b = (alf/bet) - r2
  llgrad = -c(d1a,d1b)
 } 
  
llhess  = function(phi) {
  alf = phi[1]
  bet = phi[2]
  d2aa = -trigamma(alf)
  d2ab = 1/bet
  d2bb = -alf/(bet^2)
  llhess = -matrix(c(d2aa,d2ab,d2ab,d2bb),2,2)
 } 

#GENERATE ORIGINAL DATA 
yorig1 = rgamma(n,alftru1,rate=bettru1)
yorig2 = rgamma(n,alftru2,rate=bettru2)

#COMPUTE PARAMETER ESTIMATES FOR THE ORIGINAL DATA
#GROUP 1
ycur = yorig1
ybar = mean(ycur)
yvar = var(ycur)
bet0 = ybar/yvar
alf0 = ybar*bet0
mle = optim(c(alf0,bet0),llik,gr=llgrad,lower=c(1e-4,1e-4),method="L-BFGS-B")
alfhat1 = mle$par[1]
bethat1 = mle$par[2]
exptdhat1 = alfhat1/bethat1
#GROUP 2
ycur = yorig2
ybar = mean(ycur)
yvar = var(ycur)
bet0 = ybar/yvar
alf0 = ybar*bet0
mle = optim(c(alf0,bet0),llik,gr=llgrad,lower=c(1e-4,1e-4),method="L-BFGS-B")
alfhat2 = mle$par[1]
bethat2 = mle$par[2]
exptdhat2 = alfhat2/bethat2
thhat = exptdhat1/exptdhat2
lthhat = log(thhat)

#ASYMPTOTIC CI
#FIRST FIND CI FOR LOG theta AND THEN TAKE exp
a = alfhat1
b = bethat1
hess = llhess(c(a,b))
invhes = solve(hess)
grth = c((1/a),-(1/b))
v1 = t(grth) %*% invhes %*% grth / n
a = alfhat2
b = bethat2
hess = llhess(c(a,b))
invhes = solve(hess)
grth = c((1/a),-(1/b))
v2 = t(grth) %*% invhes %*% grth / n
sddif = sqrt(v1+v2)
cihw = za*sddif
cifac = exp(cihw)
thlo = thhat/cifac
thhi = thhat*cifac
asyci = c(thlo,thhi)

#BOOTSTRAP
tstar = rep(0,M)
for (m in 1:M) {
  #GROUP 1
  ycur = rgamma(n,alfhat1,rate=bethat1)
  ybar = mean(ycur)
  yvar = var(ycur)
  bet0 = ybar/yvar
  alf0 = ybar*bet0
  mle = optim(c(alf0,bet0),llik,gr=llgrad,lower=c(1e-4,1e-4),method="L-BFGS-B")
  alfhst1 = mle$par[1]
  bethst1 = mle$par[2]
  exptdst1 = alfhst1/bethst1
  #GROUP 2
  ycur = rgamma(n,alfhat2,rate=bethat2)
  ybar = mean(ycur)
  yvar = var(ycur)
  bet0 = ybar/yvar
  alf0 = ybar*bet0
  mle = optim(c(alf0,bet0),llik,gr=llgrad,lower=c(1e-4,1e-4),method="L-BFGS-B")
  alfhst2 = mle$par[1]
  bethst2 = mle$par[2]
  exptdst2 = alfhst2/bethst2
  thhst = exptdst1/exptdst2
  lthhst = log(thhst)
  #COMPUTE ASYMPTOTIC VARIANCE FOR BOOTSTRAP SAMPLE
  a = alfhst1
  b = bethst1
  hess = llhess(c(a,b))
  invhes = solve(hess)
  grth = c((1/a),-(1/b))
  v1 = t(grth) %*% invhes %*% grth / n
  a = alfhst2
  b = bethst2
  hess = llhess(c(a,b))
  invhes = solve(hess)
  grth = c((1/a),-(1/b))
  v2 = t(grth) %*% invhes %*% grth / n
  sddif.boot = sqrt(v1+v2)
  #STUDENTIZED STATISTIC
  tstar[m] = (lthhst - lthhat)/sddif.boot
 }
q = quantile(tstar,c(0.975,0.025),type=8)
lth.lo = lthhat - q[1]*sddif
lth.hi = lthhat - q[2]*sddif
bootci = exp(c(lth.lo,lth.hi))
print(cbind(thtru,thhat))
print(cbind(asyci,bootci))

###########################################################

#SAMPLE OUTPUT

# > print(cbind(thtru,thhat))
# thtru    thhat
# [1,]     4 4.555845
# > print(cbind(asyci,bootci))
# asyci   bootci
# [1,] 2.531446 2.263915
# [2,] 8.199157 8.490438
