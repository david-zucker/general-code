# PART 2 OF R CODE FOR TARGIL 6 QUESTION 2

#READ DATA AND COMPUTE HOPT
indat1 = read.table("ex6dat.txt")
x = indat1[,1]
n = length(x)
xmin = min(x)
xmax = max(x)
xval = xmin + (xmax-xmin)*c(0:511)/511
ssdint = 1.46e-6
kerfac = 0.5/sqrt(pi)
hopt = (kerfac/(n*ssdint))^(1/5)

#NORMAL KERNEL
krnl = function(u) {
  krnl = dnorm(u)
 }

#KERNEL ESTIMATE OF DENSITY
kernest = function(h,xval) {
  nval = length(xval)
  f = rep(0,nval)
  for (i in 1:nval) {
    xv = xval[i]
    f[i] = sum(krnl((x-xv)/h))/(n*h)
   }
  return(f) 
 }   

denest1 = kernest(hopt,xval)
den=density(x, bw = "SJ", from=xmin, to=xmax, n=512)
hoptSJ = den$bw
hfac = hopt/hoptSJ
print(cbind(hopt,hoptSJ,hfac))
denest2 = den$y
dmax1 = max(denest1)
dmax2 = max(denest2)
dmax = max(dmax1,dmax2)
plot(xval,denest1,type="l",ylim=c(0,dmax))
lines(xval,denest2,lty=1,col="green")
