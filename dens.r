# SETUP
library(ks)
a = read.table('sbp1dat.txt')
xdat = a[,1]
n = length(xdat)
xmin = min(xdat)
xmax = max(xdat)
xrng = xmax - xmin
xval = xmin + xrng*c(0:500)/500

#DENSITY ESTIMATE
bw0 = hpi(xdat, deriv.order=0)
den = kdde(xdat, h=bw0, deriv.order=0, eval.points=xval)
f = den$estimate
plot(xval,f,type="l")

# 1ST DERIVATIVE ESTIMATE
bw1 = hpi(xdat, deriv.order=1)
dd1 = kdde(xdat, h=bw1, deriv.order=1, eval.points=xval)
fd1 = dd1$estimate
plot(xval,fd1,type="l")
lines(c(xmin,xmax),c(0,0),lty=2)

# 2ND DERIVATIVE ESTIMATE
bw2 = hpi(xdat, deriv.order=2)
dd2 = kdde(xdat, h=bw2, deriv.order=2, eval.points=xval)
fd2 = dd2$estimate
plot(xval,fd2,type="l")
lines(c(xmin,xmax),c(0,0),lty=2)

# ESTIMATED BIAS OF DENSITY ESTIMATE
bias = 0.5 * fd2 * (bw0^2)

# BIAS-CORRECTED DENSITY ESTIMATE
forig = f
f = f - bias

# ORIGINAL VS. CORRECTED ESTIMATE
plot(xval,f,type="l")
lines(xval,forig,lty=2)

# DENSITY ESTIMATE WITH POINTWISE CONFIDENCE BANDS
knrm = sqrt(0.5/sqrt(pi))
vihw = 1.96*knrm*sqrt(f)/sqrt(n*bw0)
fl = f - vihw
fu = f + vihw
minf = min(fl)
maxf = max(fu)
plot(xval,f,type="l",ylim=c(minf,maxf))
lines(xval,fl,lty=2)
lines(xval,fu,lty=2)
