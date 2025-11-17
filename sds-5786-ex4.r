#STATISTICS FOR DATA SCIENCE - EXERCISE SET 4

set.seed(5074283)

#QUESTION 2

fk = function(x,a,b) 1-(1-x^a)^b
fk_inv = function (u,a,b) x = (1 - (1-u)^(1/b))^(1/a)

my.edf = function(x) {
  x1 = sort(unique(x))
  n1 = length(x1)
  Fn = rep(0,n1)
  for (j in 1:n1) {
    Fn[j] = mean(x <= x1[j])
  }
  ans = list(x1=x1, Fn=Fn)
  return(ans)
}

a = 2
b = 5
n = 1000
u = runif(n)
x = fk_inv(u,a,b)
us = sort(u)
xs = sort(x)
plot(xs,us,type='l')
edf_x = my.edf(x)
lines(xs,edf_x$Fn,col='red')

###########################################################

#QUESTION 3

#theoretical distribution function
fxq2 = function(x) {
  nobs = length(x)
  ans = rep(0,nobs)
  for (i in 1:nobs) {
    xcur = x[i]
    if (xcur<0) ans[i] = 0
    if ((xcur>=0) & (xcur<=1)) ans[i] = xcur/2
    if ((xcur>=1) & (xcur<=2)) ans[i] = 1/2
    if ((xcur>=2) & (xcur<=3)) ans[i] = 1/2 + (xcur-2)/2
    if (xcur>3) ans[i] = 1
  }  
  return(ans)
}

#plot the theoretical distribution function
xtheo = seq(0,4,0.01)
fxtheo = fxq2(xtheo)
plot(xtheo,fxtheo,type='l')

#inverse distribution function
fxq2inv = function(u) {
  nobs = length(u)
  ans = rep(0,nobs)
  for (i in 1:nobs) {
    ucur = u[i]
    if ((ucur>0) & (ucur<0.5)) ans[i] = 2*ucur
    if (ucur==0.5) ans[i] = 1
    if ((ucur>0.5) & (ucur<=1)) ans[i] = 2*ucur + 1
  }
  return(ans)
}

#plot of inverse distribution function
u = seq(0.01,0.5,0.01)
q = fxq2inv(u)
plot(u,q,type='l',xlim=c(0,1),ylim=c(0,4))
points(0.5,1,pch=19)
u = c(0.501,seq(0.51,1,0.01))
q = fxq2inv(u)
lines(u,q)
points(0.5,2)

#generate sample of size 1000 from the distribution
n = 1000
u = runif(n)
x = fxq2inv(u)

#plot the theoretical distribution function
#and the empirical distribution function
plot(xtheo,fxtheo,type='l')
xs = sort(x)
edf_x = my.edf(x)
lines(xs,edf_x$Fn,col='red')

###########################################################

#QUESTION 4

#preliminaries
x = seq(0,20,0.02)
d = dgamma(x,6.5,1)
plot(x,d,type='l',ylim=c(0,0.40))
d1 = dgamma(x,8,1)
lines(x,d1,col='red')
c = 2
d2 = c*d1
lines(x,d2,col='green')
summary(d2-d)

# draw the sample
nobs = 1000
c = 2
z = NULL
for (i in 1:nobs) {
  flg = 0
  while (flg == 0) {
    u0 = runif(6)
    y = - log(1-u0)
    w = sum(y)
    u = runif(1)
    bound = dgamma(w,6.5,1) / (c*dgamma(w,6,1))
    if (u <= bound) {
      z = c(z,w)
      flg = 1
    }
  }
}

# plot the gamma density
x = c(0:2000)/100
fx = dgamma(x,6.5,1)
plot(x,fx,type='l',ylim=c(0,0.2))

# estimate the density based on the sample and plot
dens=density(z, bw='SJ')
lines(dens$x,dens$y, col='red')

#qq plot
qqrange = ((1:nobs) - 0.375)/(nobs+0.25)
qtheo = qgamma(qqrange,6.5,1)
zs = sort(z)
plot(qtheo,zs)
lines(c(0,20),c(0,20),col='red')

###########################################################

#QUESTION 5

#exact calculation
exact = ppois(10,10)
print(exact)

#simulation
nreps = 100000
xbar_vec = rep(0,nreps)
for (i in 1:nreps) {
  datcur = rpois(5,2)
  xbar_vec[i] = mean(datcur)
}
simprb = mean((xbar_vec <= 2))
print(simprb)
print(exact-simprb)


    
  