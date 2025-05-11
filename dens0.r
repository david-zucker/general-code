# UNIVARIATE KERNEL DENSITY ESTIMATION
# BANDWIDTH SELECTION BY SHEATHER AND JONES (1991) METHOD

a = read.table('sbp1dat.txt')
xdat = a[,1]
adj = 1
den = density(xdat,bw="SJ",adjust=adj)
x = den$x
f = den$y
fmin = min(f)
fmax = max(f)
histo = hist(xdat, plot=FALSE)
hden = histo$den
hmin = min(hden)
hmax = max(hden)
dlo = min(fmin,hmin)
dhi = max(fmax,hmax)
hist(xdat, xlim=range(x), ylim=c(dlo,dhi), xlab="x", ylab="density", 
 breaks='Scott', probability=TRUE)
lines(x,f,lty=1)
rug(jitter(xdat))
bw=den$bw
print(bw)
