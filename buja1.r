#Buja example

library(MASS)
library(mvtnorm)

set.seed(6543)

#SETUPS
nrep = 500
coef.fac = 3
n = 250
p = 10
xmean = 5
xsd = 0.6
corr = 0.5
p1 = p+1
df = n - p1
COEF = coef.fac/p
tstat.vec = rep(0,nrep)
pval.vec = rep(0,nrep)
keep.vec = rep(0,nrep)
kstep = log(n)

print(cbind(n,p,xmean,xsd,coef.fac))

frmla = 'Y ~ 0 + .'

zro = rep(xmean,p1)
sigmat = (1-corr)*diag(p1) + matrix(corr,p1,p1)
sigmat = (xsd^2) * sigmat

#SIMULATION LOOP

for (irep in 1:nrep) {
  
if ((irep %% 50) == 0) print(irep)

#GENERATE DATA
x = rmvnorm(n, mean=zro, sigma=sigmat)
dtfrm = as.data.frame(x)
trubet = c(0,rep(COEF,p))
dtfrm$Y = x %*% trubet + rnorm(n,0,1) 
#coefficent of x1 is 0
#all other regression coefficients are equal to COEF

#CARRY OUT VARIABLE SELECTION USING BIC

# specify min and max models 
reg.min = lm(Y ~ 0 + V1, data=dtfrm)
reg.full = lm(frmla, data=dtfrm)

#variable selection using stepAIC
reg.BIC = stepAIC(reg.min, scope=list(lower=formula(reg.min), 
  upper=formula(reg.full)), direction='both', k=kstep, trace=0) #k=kstep gives BIC

#store t statistic and p-value, and number of variables kept
tstat.vec[irep] = summary(reg.BIC)[[4]][1,3]
pval.vec[irep] = summary(reg.BIC)[[4]][1,4]
keep.vec[irep] = length(reg.BIC$coefficients)-1

}
#END OF SIMULATION LOOP

#DISTRIBUTION OF T VALUES

#HISTOGRAM
hist(tstat.vec,breaks='Scott',freq=F,main=NULL,ylim=c(0,0.4))
tt = seq(-3,3,length.out=100)
ff = dt(tt,df)
lines(tt,ff,col='red')

#QQ PLOT
prbs = (c(1:nrep) - 0.5)/nrep
qntl = qt(prbs,df)
empqntl = sort(tstat.vec)
plot(qntl,empqntl)
lines(qntl,qntl,col='red')

#COVERAGE RATE
coverage.rate = mean((pval.vec > 0.05))
print(paste('coverage rate = ', coverage.rate))

print('summary of t statistics')
print(summary(tstat.vec))
print('distribution of number of additional variables in the model')
print(table(keep.vec)/nrep)

