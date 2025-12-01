# PART 1 OF R CODE FOR TARGIL 6 QUESTION 2

#READ DATA AND COMPUTE HOPT2
indat1 = read.table("ex6dat.txt")
x = indat1[,1]
n = length(x)
xmean = mean(x)
sdx = sd(x)
iqrsd = IQR(x)/1.349
sig = min(sdx,iqrsd)
hopt2 = 0.9397*sig*(n^(-1/9))

#SECOND DERIVATIVE OF NORMAL KERNEL
krnld2 = function(u) {
  krnld2 = (u^2-1)*dnorm(u)
 }
 
#KERNEL ESTIMATE OF SECOND DERIVATIVE OF DENSITY
kernd2est = function(h,xval) {
  nval = length(xval)
  f = rep(0,nval)
  for (i in 1:nval) {
    xv = xval[i]
    f[i] = sum(krnld2((x-xv)/h))/(n*(h^3))
   }
  return(f) 
 }   

#COMPUTATION OF INTEGRAL OF SQUARED SECOND DERIVATIVE
ssdint = function(J) {
  twoJ = 2*J
  y = (2*c(1:J)-1)/twoJ
  xval = xmean + sdx*(log(y) - log(1-y))
  fd2 = kernd2est(hopt2,xval)
  g = fd2^2/(y*(1-y))
  ssdint = sdx*sum(g)/twoJ
 } 

#COMPUTATION OF THE INTEGRAL FOR VARIOUS VALUES OF J
ssdires = NULL
M = 30
jvec = 10 * c(1:M)
for (m in 1:M) {
  J = jvec[m]
  ssdi = ssdint(J)
  ssdires = c(ssdires,ssdi)
 } 

print(ssdires)  

