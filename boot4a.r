#PROGRAM TO DEMONSTRATE TWO-SAMPLE NONPARAMETRIC BOOTSTRAP
#PARAMETER OF INTEREST IS RATIO OF EXPECTATIONS

#PRELIMINARIES
set.seed(342345)
za = -qnorm(0.025)
n = 20
alftru1 = 2
bettru1 = 3
exptd1 = alftru1*bettru1
alftru2 = 3
bettru2 = 5
exptd2 = alftru2*bettru2
thtru = exptd1/exptd2
M = 1000

#GENERATE ORIGINAL DATA 
yorig1 = rgamma(n,alftru1,scale=bettru1)
yorig2 = rgamma(n,alftru2,scale=bettru2)

#COMPUTE PARAMETER ESTIMATE FOR THE ORIGINAL DATA
theta_hat = mean(yorig1)/mean(yorig2)

#BOOTSTRAP
tstar = rep(0,M)
for (m in 1:M) {
  ix1 = sample.int(n,n,replace=TRUE)
  exp1star = mean(yorig1[ix1])
  ix2 = sample.int(n,n,replace=TRUE)
  exp2star = mean(yorig2[ix2])
  theta_hat_star = exp1star/exp2star
  tstar[m] = theta_hat_star - theta_hat
}  
q = quantile(tstar,c(0.975,0.025))
bootci = theta_hat*exp(-q)
print(cbind(thtru,theta_hat))
print(bootci)  

