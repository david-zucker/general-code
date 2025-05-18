# 2D KERNEL DENSITY ESTIMATION
# BANDWIDTH SELECTION USING METHOD OF WAND AND JONES (1994)

# SETUP
library(ks)
library(lattice)
npt = 80
np1 = npt + 1
np2 =  np1^2
a = read.table('birthdat.txt')
x = a[,1]
y = a[,2]
n = length(x)
datmat = cbind(x,y)
xmin = 0.75*min(x)
xmax = 1.25*max(x)
xmean = mean(x)
xrng = xmax - xmin
xval = xmin + xrng*c(0:npt)/npt
ymin = 0.75*min(y)
ymax = 1.25*max(y)
ymean = mean(y)
yrng = ymax - ymin
yval = ymin + yrng*c(0:npt)/npt
eval = NULL
for (i in 1:np1) {
  xv = rep(xval[i],np1)
  ev = cbind(xv,yval)
  eval = rbind(eval,ev)
 }

# NORMALIZE DATA
sdx = sd(x)
sdy = sd(y)
nrmfac1 = cbind(rep(sdx,n),rep(sdy,n))
nrmfac2 = cbind(rep(sdx,np2),rep(sdy,np2))
dm2 = datmat/nrmfac1
ev2 = eval/nrmfac2

# UNIVARIATE DENSITY ESTIMATES FOR X AND Y
# (NORMALIZED AND CENTERED)
xdnew = (x-xmean)/sdx
dx = density(xdnew,bw="SJ")
plot(dx$x,dx$y,type="l")
ydnew = (y-ymean)/sdy
dy = density(ydnew,bw="SJ")
plot(dy$x,dy$y,type="l")

# DENSITY ESTIMATE
Hmat = Hpi(dm2, deriv.order=0)
den = kdde(dm2, H=Hmat, deriv.order=0, eval.points=ev2)
f0 = den$estimate
f = den$estimate/(sdx*sdy)
f = 1e5 * f

# INFO ON H MATRIX
egn = eigen(Hmat)
evals = egn$values
evec = egn$vectors
a = sqrt(evals)
arat = a[2]/a[1]
theta = atan(evec[2,1]/evec[1,1])
if (theta < 0) {theta = theta + 2*pi}
theta = 360*theta/(2*pi)
Hmat
evals
evec
a
arat
theta

# LEVEL PLOT OF DENSITY ON NORMALIZED SCALE
library(graphics)
zv0 = matrix(0,np1,np1)
for (k1 in 1:np1) {
for (k2 in 1:np1) {
  ix = np1*(k1-1) + k2
  zv0[k1,k2] = f0[ix]
 }
 }
lvl = seq(from=0, to=0.2, by=0.02)
clrs = c("white","yellow","lightgreen","lavender","darkblue",
  "salmon","darkgreen","darkblue","purple","darkred")
shft = -0.842  
filled.contour(x=(xval-xmean)/sdx, y=(yval-ymean)/sdy, z=zv0, xlim = c(-3,3),col=clrs,
  ylim = c(-3,3), zlim = c(0,0.2), levels=lvl)
lines(shft+5*c(0,evec[1,1]),5*c(0,evec[2,1]),lty=1,col="black")
lines(shft-5*c(0,evec[1,1]),-5*c(0,evec[2,1]),lty=1,col="black")
lines(shft+5*c(0,evec[1,2]),5*c(0,evec[2,2]),lty=1,col="black")
lines(shft-2.5*c(0,evec[1,2]),-2.5*c(0,evec[2,2]),lty=1,col="black")
lines(c(shft,shft),c(-3,3),lty=2)
lines(c(-3.2,1.5),c(0,0),lty=2)

# LEVEL PLOT OF DENSITY ON ORIGINAL SCALE
library(graphics)
zv = matrix(0,np1,np1)
for (k1 in 1:np1) {
for (k2 in 1:np1) {
  ix = np1*(k1-1) + k2
  zv[k1,k2] = f[ix]
 }
 }
lvl = seq(from=0, to=2.25, by=0.25)
clrs = c("white","yellow","lightgreen","lavender","darkblue",
  "salmon","darkgreen","darkblue","purple","darkred")
filled.contour(x=xval, y=yval, z=zv, xlim = c(20,110),col=clrs,
  ylim = c(500,5000), zlim = c(0,2.5), levels=lvl)

# PREPARE FOR WIREFRAME PLOTTING
xv = eval[,1]
yv = eval[,2]
df = data.frame(xv,yv,eval)
scl = list(arrows=FALSE, tick.number = 5)

# WIREFRAME PLOTTING FUNCTION
denplt = function(x_ang, y_ang) {
  scr = list(x=x_ang, y=y_ang)
  wireframe(f~xv*yv,data=df,drape=TRUE, xlab = "MW", ylab = "IW", zlab = "f",
    main = "2D Kernel Density Estimate", xlim = range(xval), ylim = range(yval),
    zlim = range(f), scales=scl, screen=scr)
  }  
  
# DENSITY SURFACE PLOTS  
denplt(-30,10)  
denplt(-40,0)


###############################################################################

# OUTPUT FROM INFO ON H MATRIX
# 
# > Hmat
#            [,1]       [,2]
# [1,] 0.09316098 0.02951363
# [2,] 0.02951363 0.16762598
# > evals
# [1] 0.17790468 0.08288229
# > evec
#           [,1]       [,2]
# [1,] 0.3288940 -0.9443668
# [2,] 0.9443668  0.3288940
# > a
# [1] 0.4217875 0.2878929
# > arat
# [1] 0.6825543
# > theta
# [1] 70.79834
