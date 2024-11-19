


fxn = function(x,n) {     
  kvec = 0:floor(x)
  ktrm = (-1)^kvec * choose(n,kvec) * (x-kvec)^n 
  ans = sum(ktrm)/factorial(n)
  return(ans)
}  

nvec = c(10, 15, 20, 25)
epsvec = c(0.01, 0.05, 0.10, 0.15, 0.20)

for (ixn in 1:4) {
  
  n = nvec[ixn]
  print.noquote(paste0('n = ', n))
  ansr = NULL
  
  for (ixe in 1:3) {
    eps = epsvec[ixe]
    eps2 = eps^2
    pi0 = 1 - fxn(0.5*n*(eps+1),n) + fxn(0.5*n*(-eps+1),n)
    pi1 = 1/(3*n*eps2)
    pi2.arg = 1/sqrt(pi1)
    pi2 = 1 - pnorm(pi2.arg) + pnorm(-pi2.arg)
    pi3 = 2 * exp(-0.5*eps2*n)
    anscur = c(pi0,pi1,pi2,pi3)
    ansr = rbind(ansr,anscur)
  }  

  ansr = as.data.frame(ansr)
  colnames(ansr) = c('pi0','pi1','pi2','pi3')
  rownames(ansr) = c('eps=0.01','eps=0.05', 'eps=0.10')
  print(ansr)
  cat('\n')
  
}

# > source("C:/Users/owner/Desktop/ex3.r")
# [1] n = 10
# pi0        pi1       pi2      pi3
# eps=0.01 0.9569785 333.333333 0.9563199 1.999000
# eps=0.05 0.7872974  13.333333 0.7841912 1.975156
# eps=0.10 0.5890374   3.333333 0.5838824 1.902459
# 
# [1] n = 15
# pi0        pi1       pi2      pi3
# eps=0.01 0.9470523 222.222222 0.9465164 1.998501
# eps=0.05 0.7397674   8.888889 0.7373157 1.962849
# eps=0.10 0.5060142   2.222222 0.5023350 1.855487
# 
# [1] n = 20
# pi0        pi1       pi2      pi3
# eps=0.01 0.9387208 166.666667 0.9382579 1.998001
# eps=0.05 0.7005897   6.666667 0.6985354 1.950620
# eps=0.10 0.4413567   1.666667 0.4385780 1.809675
# 
# [1] n = 25
# pi0        pi1       pi2      pi3
# eps=0.01 0.9314007 133.333333 0.9309874 1.997502
# eps=0.05 0.6667839   5.333333 0.6650055 1.938466
# eps=0.10 0.3886361   1.333333 0.3864762 1.764994

    
    
  
