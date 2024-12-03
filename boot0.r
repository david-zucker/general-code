#DEMO OF USING SIMULATION FOR PROBABILITY CALCULATION

#gamma distribution with Wasserman parameterization

#SETUP
set.seed(342345)
n = 20
alftru = 0.25
bettru = 1
M = 100000
c = 0.225

#EXACT CALCULATION
alfnew = n*alftru
betnew = bettru/n
exactp = pgamma(c,alfnew,scale=betnew)

#ASYMPTOTIC CALCULATION
expval = alftru*bettru
varval = alftru*(bettru^2)
z = sqrt(n)*(c-expval)/sqrt(varval)
asyp = pnorm(z)

#SIMULATION CALCULATION
ybar = rep(0,M)
for (m in 1:M) {
  y = rgamma(n,alftru,scale=bettru)
  ybar[m] = mean(y)
 }
indic = (ybar <= c)
psim = sum(indic)/M

#OUTPUT
print(cbind(exactp,asyp,psim))

###########################################################

#SAMPLE OUTPUT

#        exactp      asyp    psim
#[1,] 0.4678964 0.4115316 0.46794
