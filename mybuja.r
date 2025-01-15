#modified Buja example
#selection based on p-values of the coefficients instead of BIC

set.seed(6543)

#SETUPS
nrep = 5000
coef.fac = 0
n = 250
p = 10
xmean = 1
xsd = 0.6
pcut = 0.15
p1 = p+1
df = n - p1
COEF = coef.fac/p
tstat.vec = NULL
keep.vec = NULL

print(cbind(n,p,xmean,xsd,coef.fac,pcut))

for (irep in 1:nrep) {

#GENERATE DATA
x = rnorm(p1*n,xmean,xsd)
x = matrix(x,n,p1)
trubet = c(0,rep(COEF,p))
y = x %*% trubet + rnorm(n,0,1) 
#coefficent of x1 is 0
#all other regression coefficients are equal to COEF

#PERFORM REGRESSION INCLUDING ALL VARIABLES
xtx = t(x) %*% x
xtx.i = solve(xtx)
xty = t(x) %*% y
beta = xtx.i %*% xty
prd = x %*% beta
res = y - prd
sig2 = sum(res^2)/df
vars = sig2*diag(xtx.i)
sds = sqrt(vars)
tstat = beta/sds
pval = 2*(1 - pt(abs(tstat),df))

#SELECT VARIABLES WITH P-VALUE <= pcut
pval1 = pval[2:p1]
keep = 1 + which(pval1 <= pcut)
keep.vec = c(keep.vec,length(keep))
keep = c(1,keep)

#RERUN REGRESSION
xnew = x[,keep]
xtx.new = t(xnew) %*% xnew
xtx.i.new = solve(xtx.new)
xty.new = t(xnew) %*% y
beta.new = xtx.i.new %*% xty.new
prd.new = xnew %*% beta.new
res.new = y - prd.new
sig2.new = sum(res.new^2)/df

#COMPUTE T STATISTIC
tstat.cur = beta.new[1]/sqrt(sig2.new*xtx.i.new[1,1])
tstat.vec = c(tstat.vec,tstat.cur)

}

hist(tstat.vec,breaks='Scott',freq=F,main=NULL,ylim=c(0,0.4))
tt = seq(-3,3,length.out=100)
ff = dt(tt,df)
lines(tt,ff,col='red')

tcrit = qt(0.975,df)
coverage.rate = mean((abs(tstat.vec) <= tcrit))
print(paste('coverage rate = ', coverage.rate))

print('summary of t statistics')
print(summary(tstat.vec))
print('distribution of number of additional variables in the model')
print(table(keep.vec)/nrep)

#        n  p xmean xsd coef.fac pcut
# [1,] 250 10     1 0.6        0 0.15
# [1] "coverage rate =  0.8894"
# [1] "summary of t statistics"
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -3.932692 -0.891387 -0.005587 -0.005580  0.867558  3.711393 
# [1] "distribution of number of additional variables in the model"
# keep.vec
#      0      1      2      3      4      5      6      7 
# 0.2132 0.3336 0.2674 0.1328 0.0408 0.0096 0.0022 0.0004 
# > 

#        n  p xmean xsd coef.fac pcut
# [1,] 250 10     1 0.6    -0.35 0.15
# [1] "coverage rate =  0.761"
# [1] "summary of t statistics"
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -7.6212 -1.7715 -0.7138 -1.0412  0.1714  3.5050 
# [1] "distribution of number of additional variables in the model"
# keep.vec
#      0      1      2      3      4      5      6      7 
# 0.1340 0.3338 0.3030 0.1574 0.0558 0.0134 0.0024 0.0002 
