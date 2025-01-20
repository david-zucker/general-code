indat = read.table("fatdat.txt")
names(indat) = c('id', 'fatbro', 'fatsir', 'density', 'age', 'wt', 'ht', 
  'BMI', 'ffw', 'neck', 'chest', 'abdo','hip', 'thigh', 'knee', 'ankle', 
  'biceps', 'forearm', 'wrist')

#the y variable is fatbro: percent body fat using Brozek's method
#we will use the following x variables: age, adip, chest, abdo, hip, thigh)

#BMI = wt.kg/(ht.m^2)

#PART (a)
reg = lm(fatbro ~ age + abdo, data=indat)
print(summary(reg))

#PART (b)
n = dim(indat)[1]
wuns = rep(1,n)
x = cbind(wuns,indat$age,indat$abdo)
betvec = solve(t(x)%*%x) %*% t(x) %*% indat$fatbro
print(betvec)
cat('\n')

#PART (c)
tstat = 0.06605 / 0.02290
print(tstat)
pval = 2*(1 - pt(tstat,249))
print(pval)
cat('\n')

#PART (d)
tcrit = qt(0.975,249)
print(tcrit)
cat('\n')
cihw = tcrit*0.0268
b_lo = 0.5671 - cihw
b_hi = 0.5671 + cihw
print(cbind(b_lo,b_hi))
cat('\n')

#PART (e)
a = c(1,30,80)
bet_est = coef(reg)
th_est = sum(a*bet_est)
print(th_est)
cat('\n')
covmat = vcov(reg)
print(covmat)
varth = t(a) %*% covmat %*% a
print(varth)
cihw = tcrit * sqrt(varth)
th_lo = th_est - cihw
th_hi = th_est + cihw
print(cbind(th_lo,th_hi))
cat('\n')

#PART (f)
ee = residuals(reg)
px = x %*% solve(t(x)%*%x) %*% t(x)
pxdiag = diag(px)
mse_cv = mean((ee/(1-pxdiag))^2)
print(mse_cv)



 