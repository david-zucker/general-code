#PROGRAM TO DEMONSTRATE NONPARAMETRIC DOUBLE STUDENTIZED BOOTSTRAP
#PARAMETER OF INTEREST IS COEFFICIENT OF VARIATION E[Y]/SD(Y)

#PRELIMINARIES
set.seed(342345)
n = 20
alfgamtru = 3
betgamtru = 0.2
thtru = sqrt(alfgamtru)
R = 1000
M = 1000
alf = 0.05
nomcov = 1-alf
za = -qnorm(alf/2)
input = cbind(n,R,M,thtru,alf)

#DELTA METHOD SD OF theta-hat
delsd = function(ycur) {
  wcur = ycur^2
  ucur = cbind(ycur,wcur)
  ucov = cov(ucur)
  u1 = mean(ycur)
  u2 = mean(wcur)
  fac1 = u2-u1^2
  fac2 = fac1^1.5
  grth1 = (fac1 + u1^2)/fac2
  grth2 = -0.5*u1/fac2
  grth = c(grth1,grth2)
  delvar = t(grth) %*% ucov %*% grth / n
  delsd = sqrt(delvar)
 } 

#BOOTSTRAP ESTIMATE OF COVERAGE PROBABILITIES
covest = function(eta) {
  covest = 0
  for (r in 1:R) {
    wststcur = wstst[r,]
    plo = eta/2
    phi = 1-plo
    q = quantile(wststcur,c(plo,phi),type=8)
    q = as.vector(q)
    cilo = thstar[r] - q[2]*sdthb[r]
    cihi = thstar[r] - q[1]*sdthb[r]
    covest = covest + (cilo <= thhat)*(thhat <= cihi)
   } 
  covest = covest/R
  return(covest)
 }   
    
#OBJECTIVE FUNCTION FOR eta SEARCH
objfn = function(eta) {
  covr = covest(eta)
  objfn = (covr - nomcov)^2
 } 

#GENERATE SAMPLE OF SIZE n FROM GAMMA(alfgamtru,betgamtru) DISTRIBUTION
yorig = rgamma(n,alfgamtru,scale=betgamtru)

#COMPUTE ESTIMATE OF THETA FROM THE ORIGINAL DATA
ycur = yorig
thhat = mean(ycur)/sd(ycur)
sdth = delsd(ycur)

#ASYMPTOTIC CI
acihw = za*sdth
acilo = thhat - acihw
acihi = thhat + acihw

#BOOTSTRAP QUANTITIES
thstar = rep(0,R)
wstar = rep(0,R)
sdthb = rep(0,R)
wstst = matrix(0,R,M)
for (r in 1:R) {
  yboot = sample(yorig, size=n, replace=TRUE) 
  ycur = yboot
  thhst = mean(ycur)/sd(ycur)
  sdthboot = delsd(ycur)
  thstar[r] = thhst
  sdthb[r] = sdthboot
  wstar[r] = (thhst - thhat)/sdthboot
  for (m in 1:M) { 
    ycur = sample(yboot, size=n, replace=TRUE) 
    thhstst = mean(ycur)/sd(ycur)
    sdthcur = delsd(ycur)
    wstst[r,m] = (thhstst - thhst)/sdthcur
   } 
 }
 
#SEARCH FOR eta
etsr = optimize(objfn,c(1e-5,1-1e-5))
etaopt = etsr$minimum

#COMPUTATION OF FINAL CONFIDENCE INTERVAL 
prlo = etaopt/2
prhi = 1-prlo
q = quantile(wstar,c(prlo,prhi),type=4)
q = as.vector(q)
cilo = thhat - q[2]*sdth
cihi = thhat - q[1]*sdth
print(input)
print(cbind(thhat,sdth,acilo,acihi,etaopt,q[1],q[2],cilo,cihi)) 
 
###########################################################

#SAMPLE OUTPUT (SLIGHTLY EDITED)

#  n    R    M    thtru  alf
# 20 1000 1000 1.732051 0.05
#
#    thhat     sdth     acilo    acihi     etaopt   
# 2.285347 0.2682499 1.759587 2.811107 0.05192002 
#
#       q[1]     q[2]     cilo     cihi 
#  -1.977693 1.869713 1.783797 2.815863


#   -1.809 1.7601220 1.813195 2.770611
