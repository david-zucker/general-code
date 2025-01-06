#PROGRAM TO COMPUTE MAXIMUM LIKELIHOOD ESTIMATES FOR THE WEIBULL DISTRIBUTION
#AND POINTWISE CONFIDENCE INTERVAL FOR THE CUMULATIVE DISTRIBUTION FUNCTION
#INCLUDES ALSO MODEL CHECKING

set.seed(5074283)
cat('\n')

#READ DATA
indat = read.table("ex9dat.txt")
z = indat[,1]
n = length(z)

#WEIBULL CDF
wbcdf = function(theta,xv) {
  th1 = theta[1]
  th2 = theta[2]
  eth1 = exp(th1)
  eth2 = exp(th2)
  cdfval= 1 - exp(-(eth1*xv)^eth2)
  return(cdfval)
}  

#INVERSE WEIBULL CDF
#WEIBULL CDF
wbcdf_inv = function(theta,uu) {
  th1 = theta[1]
  th2 = theta[2]
  eth2 = exp(th2)
  invcdf = exp(log(-log(1-uu))/eth2 - th1)
  return(invcdf)
}  

#GRADIENT OF WEIBULL CDF
wbcdfgrd = function(theta,x) {
  th1 = theta[1]
  th2 = theta[2]
  eth1 = exp(th1)
  eth2 = exp(th2)
  expexp = eth1^eth2
  expexpxx = (eth1*x)^eth2
  srv = exp(-expexpxx)
  G1 = srv*eth2*expexpxx
  G2 = srv*eth2*expexpxx*(log(x)+th1)
  return(c(G1,G2))
}  

#WEIBULL DENSITY
wbdens = function(theta,x) {
  th1 = theta[1]
  th2 = theta[2]
  eth1 = exp(th1)
  eth2 = exp(th2)
  expexpxx1 = (eth1*x)^(eth2-1)
  expexpxx2 = (eth1*x)^eth2
  dens = eth1*eth2*expexpxx1*exp(-expexpxx2)
  return(dens)
}  
  
#LOG-LIKELIHOOD FUNCTION
llfun = function(theta,x) {
  dens = wbdens(theta,x)
  llval = sum(log(dens))
  return(llval)
}

#GRADIENT OF LOG-LIKELIHOOD FUNCTION
llgrd = function(theta,x) {
  th1 = theta[1]
  th2 = theta[2]
  eth1 = exp(th1)
  eth2 = exp(th2)
  expexp = eth1^eth2
  n = length(x)
  xeth2 = x^eth2 
  slogx = sum(log(x))
  sxeth2 = sum(xeth2)
  slxeth2 = sum(log(x)*xeth2)
  G1 = n*eth2 - eth2*expexp*sxeth2
  G2 = n*th1*eth2 + n + eth2*slogx - th1*eth2*expexp*sxeth2 - eth2*expexp*slxeth2 
  grdval = c(G1,G2)
  return(grdval)
}

#HESSIAN OF LOG-LIKELIHOOD FUNCTION
llhess = function(theta,x) {
  th1 = theta[1]
  th2 = theta[2]
  eth1 = exp(th1)
  eth2 = exp(th2)
  expexp = eth1^eth2
  n = length(x)
  xeth2 = x^eth2 
  slogx = sum(log(x))
  sxeth2 = sum(xeth2)
  slxeth2 = sum(log(x)*xeth2)
  sl2xeth2 = sum((log(x)^2)*xeth2)
  hess = matrix(0,2,2)
  hess[1,1] = -(eth2^2)*expexp*sxeth2 
  hess[1,2] = n*eth2 - eth2*expexp*sxeth2 - th1*(eth2^2)*expexp*sxeth2 -
    (eth2^2)*expexp*slxeth2
  hess[2,1] = hess[1,2]
  hess[2,2] = n*th1*eth2 + eth2*slogx - th1*eth2*expexp*sxeth2 - (th1^2)*(eth2^2)*expexp*sxeth2 -
    th1*(eth2^2)*expexp*slxeth2 - eth2*expexp*slxeth2 - th1*(eth2^2)*expexp*slxeth2 -
    (eth2^2)*expexp*sl2xeth2
  return(hess)
}

#SEARCH FOR STARTING VALUES
#WE HAD TO FIDDLE WITH THE SEARCH METHOD (trying different values of K and sig) 
#TO GET A GOOD STARTING VALUE
K = 50
sig = 1
cvals = rnorm(K,0,sig)
dvals = rnorm(K,0,sig)
lval = rep(0,K)
for (k in 1:K) {
  theta = c(cvals[k],dvals[k])
  lval[k] = llfun(theta,z)
}
k_opt = which.max(lval)
thtcur = c(cvals[k_opt],dvals[k_opt])

#SET UP FOR NEWTON-RAPHSON ITERATIONS
maxiter = 50
tol = 1e-6
diff = 99
iter = 0
rescur = c(0,thtcur)
rslts_mle = NULL

#CARRY OUT NEWTON-RAPHSON ITERATIONS
while ((iter < maxiter) & (diff > tol)) {
  grdval = llgrd(thtcur,z)
  rescur = c(rescur,grdval)
  rslts_mle = rbind(rslts_mle,rescur)
  iter = iter + 1
  hess = llhess(thtcur,z)
  step = solve(hess,grdval)
  thtcur = thtcur - step
  rescur = c(iter,thtcur)
  diff = max(abs(step))
}
grdval = llgrd(thtcur,z)
rescur = c(rescur,grdval)
rslts_mle = rbind(rslts_mle,rescur)
thtfin = thtcur

#PRINT RESULTS
rslts = as.data.frame(rslts_mle)
rownames(rslts_mle) = NULL
colnames(rslts_mle) = c('iter', 'th1', 'th2', 'g1', 'g2')
print(rslts_mle)
cat('\n')

#COVARIANCE MATRIX OF MLE = INVERSE INFORMATION MATRIX
vmat = -solve(llhess(thtfin,z)/n)

#CONFIDENCE INTERVAL FOR CDF AT THE POINT a=12
a=12
alf = 0.05
za = qnorm(1-alf/2)
Fhat = wbcdf(thtfin,a)
Fgrdcur = wbcdfgrd(thtfin,a)
gam2mod = t(Fgrdcur) %*% vmat %*% Fgrdcur
cihw = za*sqrt(gam2mod)/sqrt(n)
F_lo = Fhat - cihw
F_hi = Fhat + cihw
print(cbind(a,Fhat,F_lo,F_hi))

#ALTERNATE CONFIDENCE INTERVAL
th1 = thtfin[1]
th2 = thtfin[2]
eth2 = exp(th2)
alog = log(a)
arg = eth2*(th1+alog)
arggrd = rep(0,2)
arggrd[1] = eth2
arggrd[2] = arg
gam_alt = t(arggrd) %*% vmat %*% arggrd
cihw1 = za*sqrt(gam_alt)/sqrt(n)
arg_lo = arg - cihw1
arg_hi = arg + cihw1
F_lo_alt = 1 - exp(-exp(arg_lo))
F_hi_alt = 1 - exp(-exp(arg_hi))
print(cbind(a,Fhat,F_lo_alt,F_hi_alt))

#DENSITY PLOT
dens = density(z,bw='SJ')
y0 = wbdens(thtfin,dens$x)
plot(dens$x,y0,type='l')
lines(dens$x,dens$y,col='red')

#QQ PLOT
uvec = ((1:n)-0.375)/(n+0.25)
quant = wbcdf_inv(thtfin,uvec)
zsort = sort(z)
plmax = 1.1*max(c(zsort,quant))
plot(quant,zsort,ylim=c(0,plmax),type='l')
lines(c(0,plmax),c(0,plmax),col='red')


