setwd("C:/Users/owner/Dropbox/WORK/Zucker/Stat for DS/Targilim/5786")

set.seed(2295378)

#READ DATA
indat = read.table('bmi.txt')
x = indat[,1]

#remove the extreme outlier for purposes of this demonstration
ix = which(x>40)
x = x[-ix]

#NORMAL-THEORY CI WITH DELTA METHOD VARIANCE (FROM EX 3)
skci1 = function(x,clvl) {
  
  #CRITICAL VALUE
  alfhlf = (1-clvl)/2
  zcrit = qnorm(1-alfhlf)

  #COMPUTE REQUIRED SAMPLE MOMENTS
  n = length(x)
  mm = rep(0,6)
  for (j in 1:6) {
    mm[j] = mean(x^j)
  }

  #COVARIANCE MATRIX TERMS
  A = mm[2] - mm[1]^2
  B = mm[3] - mm[1]*mm[2]
  C = mm[4] - mm[2]^2
  D = mm[4] - mm[1]*mm[3]
  E = mm[5] - mm[2]*mm[3]
  F = mm[6] - mm[3]^2
  sig = sqrt(A)
  sig3 = sig^3
  sig5 = sig^5
  mu3hat = mm[3] - 3*mm[2]*mm[1] + 2*(mm[1]^3)
  gamhat = mu3hat/sig3

  #GRADIENT TERMS
  G = 3*mm[1]*mu3hat/sig5 + (6*(mm[1]^2) - 3*mm[2])/sig3
  H = -(3/2)*mu3hat/sig5 - 3*mm[1]/sig3
  K = 1/sig3
  
  #WRAP-UP
  tau2 = A*(G^2) + 2*B*G*H + 2*D*G*K + C*(H^2) + 2*E*H*K + F*(K^2)
  cihw = zcrit * (1/sqrt(n)) * sqrt(tau2)
  gam_lo = gamhat - cihw
  gam_hi = gamhat + cihw
  ans = cbind(gamhat,gam_lo,gam_hi)
  return(ans)
  
} 

#NORMAL-THEORY CI WITH BOOTSTRAP VARIANCE
skci2 = function(x,clvl,nboot) {
  
  #CRITICAL VALUE
  alfhlf = (1-clvl)/2
  zcrit = qnorm(1-alfhlf)
  
  #COMPUTE REQUIRED SAMPLE MOMENTS AND ESTIMATE ON ORIGINAL DATA
  n = length(x)
  mm = rep(0,3)
  for (j in 1:3) {
    mm[j] = mean(x^j)
  }
  sig = sqrt(mm[2]-mm[1]^2)
  mu3hat = mm[3] - 3*mm[2]*mm[1] + 2*(mm[1]^3)
  gamhat = mu3hat/sig^3
  
  #BOOTSTRAP LOOP
  tstar = rep(0,nboot)
  for (b in 1:nboot) {
    #draw a sample from the data
    xcur = sample(x,n,replace=TRUE)
    #compute estimate on bootstrap sample
    mmcur = rep(0,3)
    for (j in 1:3) {
      mmcur[j] = mean(xcur^j)
    }
    sigcur = sqrt(mmcur[2]-mmcur[1]^2)
    mu3hatcur = mmcur[3] - 3*mmcur[2]*mmcur[1] + 2*(mmcur[1]^3)
    gamhatcur = mu3hatcur/sigcur^3
    #store tstar value for current bootstrap replication
    tstar[b] = gamhatcur - gamhat
  }
  
  #WRAP-UP
  cihw = zcrit * sd(tstar)
  gam_lo = gamhat - cihw
  gam_hi = gamhat + cihw
  ans = cbind(gamhat,gam_lo,gam_hi)
  return(ans)
  
}  
  
#FULLY NONPARAMETRIC BOOTSTRAP CI
skci3 = function(x,clvl,nboot) {

  #SET UP FOR CRITICAL VALUES
  alfhlf = (1-clvl)/2
  prb1 = alfhlf
  prb2 = 1-alfhlf

  #COMPUTE REQUIRED SAMPLE MOMENTS AND ESTIMATE ON ORIGINAL DATA
  n = length(x)
  mm = rep(0,3)
  for (j in 1:3) {
    mm[j] = mean(x^j)
  }
  sig = sqrt(mm[2]-mm[1]^2)
  mu3hat = mm[3] - 3*mm[2]*mm[1] + 2*(mm[1]^3)
  gamhat = mu3hat/sig^3

  #BOOTSTRAP LOOP
  tstar = rep(0,nboot)
  for (b in 1:nboot) {
    #draw a sample from the data
    xcur = sample(x,n,replace=TRUE)
    #compute estimate on bootstrap sample
    mmcur = rep(0,3)
    for (j in 1:3) {
      mmcur[j] = mean(xcur^j)
    }
    sigcur = sqrt(mmcur[2]-mmcur[1]^2)
    mu3hatcur = mmcur[3] - 3*mmcur[2]*mmcur[1] + 2*(mmcur[1]^3)
    gamhatcur = mu3hatcur/sigcur^3
    #store tstar value for current bootstrap replication
    tstar[b] = gamhatcur - gamhat
  }

  #WRAP-UP
  q = quantile(tstar,type=4,probs=c(prb1,prb2))
  gam_lo = gamhat - q[2]
  gam_hi = gamhat - q[1]
  ans = cbind(gamhat,gam_lo,gam_hi)
  
  #LET'S TAKE A LOOK AT THE DISTRIBUTION OF THE BOOTSTAP TSTAR VALUES
  hist(tstar,breaks='Scott',freq=FALSE, main='histogram of tstar - regular bootstrap')
  
  return(ans)

}  

#STUDENTIZED BOOTSTRAP CI
skci4 = function(x,clvl,nboot) {
  
  #SET UP FOR CRITICAL VALUES
  alfhlf = (1-clvl)/2
  prb1 = alfhlf
  prb2 = 1-alfhlf
  
  #COMPUTE REQUIRED SAMPLE MOMENTS AND ESTIMATE ON ORIGINAL DATA
  n = length(x)
  mm = rep(0,6)
  for (j in 1:6) {
    mm[j] = mean(x^j)
  }
  A = mm[2] - mm[1]^2
  B = mm[3] - mm[1]*mm[2]
  C = mm[4] - mm[2]^2
  D = mm[4] - mm[1]*mm[3]
  E = mm[5] - mm[2]*mm[3]
  F = mm[6] - mm[3]^2
  sig = sqrt(A)
  sig3 = sig^3
  sig5 = sig^5
  mu3hat = mm[3] - 3*mm[2]*mm[1] + 2*(mm[1]^3)
  gamhat = mu3hat/sig3
  G = 3*mm[1]*mu3hat/sig5 + (6*(mm[1]^2) - 3*mm[2])/sig3
  H = -(3/2)*mu3hat/sig5 - 3*mm[1]/sig3
  K = 1/sig3
  tau2 = A*(G^2) + 2*B*G*H + 2*D*G*K + C*(H^2) + 2*E*H*K + F*(K^2)
  
  #BOOTSTRAP LOOP
  tstar = rep(0,nboot)
  tau2star_vec = rep(0,nboot)
  for (b in 1:nboot) {
    #draw a sample from the data
    xcur = sample(x,n,replace=TRUE)
    #compute estimate and asymptotic variance on bootstrap sample
    mmcur = rep(0,6)
    for (j in 1:6) {
      mmcur[j] = mean(xcur^j)
    }
    Acur = mmcur[2] - mmcur[1]^2
    Bcur = mmcur[3] - mmcur[1]*mmcur[2]
    Ccur = mmcur[4] - mmcur[2]^2
    Dcur = mmcur[4] - mmcur[1]*mmcur[3]
    Ecur = mmcur[5] - mmcur[2]*mmcur[3]
    Fcur = mmcur[6] - mmcur[3]^2
    sigcur = sqrt(Acur)
    sig3cur = sigcur^3
    sig5cur = sigcur^5
    mu3hatcur = mmcur[3] - 3*mmcur[2]*mmcur[1] + 2*(mmcur[1]^3)
    gamhatcur = mu3hatcur/sig3cur
    Gcur = 3*mmcur[1]*mu3hatcur/sig5cur + (6*(mmcur[1]^2) - 3*mmcur[2])/sig3cur
    Hcur = -(3/2)*mu3hatcur/sig5cur - 3*mmcur[1]/sig3cur
    Kcur = 1/sig3cur
    tau2cur = Acur*(Gcur^2) + 2*Bcur*Gcur*Hcur + 2*Dcur*Gcur*Kcur + Ccur*(Hcur^2) +
      2*Ecur*Hcur*Kcur + Fcur*(Kcur^2)
    #store tstar value for current bootstrap replication
    tstar[b] = (gamhatcur - gamhat)/sqrt(tau2cur)
    tau2star_vec[b] = tau2cur
  }
  
  #WRAP-UP
  q = quantile(tstar,type=4,probs=c(prb1,prb2))
  sd_est = sqrt(tau2)
  gam_lo = gamhat - q[2]*sd_est
  gam_hi = gamhat - q[1]*sd_est
  ans = cbind(gamhat,gam_lo,gam_hi)

  #LET'S TAKE A LOOK AT THE DISTRIBUTION OF THE BOOTSTAP TSTAR VALUES
  hist(tstar,breaks='Scott',freq=FALSE, main='histogram of tstar - studentized bootstrap')
  
  #LET'S TAKE A LOOK AT THE DISTRIBUTION OF THE BOOTSTAP TAU2 VALUES
  hist(tau2star_vec,breaks='Scott',freq=FALSE, main='histogram of tau2 values')
  plot(ecdf(tau2star_vec))
  
  return(ans)
  
}  

ci1 = skci1(x,0.95)
ci2 = skci2(x,0.95,500)
ci3 = skci3(x,0.95,10000)
ci4 = skci4(x,0.95,10000)
cat('\n')
print(noquote('Normal-theory CI with delta method variance'))
print(ci1)
cat('\n')
print(noquote('Normal-theory CI with bootstrap variance'))
print(ci2)
cat('\n')
print(noquote('Fully nonparametric bootstrap CI'))
print(ci3)
cat('\n')
print(noquote('Studentized bootstrap CI'))
print(ci4)
