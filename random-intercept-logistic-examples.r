#TOENAIL DATA
library(GLMMadaptive)
dat1 = read.table('toenaildat.txt',header=T)
print(dat1[1:10,])
dat1$trt.time = dat1$treatn * dat1$time
res1 = mixed_model(y ~ time + treatn + trt.time, ~ 1|idnew, data=dat1, family=binomial(),
  control=list(nAGQ=20))
print(summary(res1))
Omega1 = solve(res1$Hessian)
print(Omega1)

# Call:
#   mixed_model(fixed = y ~ time + treatn + trt.time, random = ~1 | 
#                 idnew, data = dat1, family = binomial(), control = list(nAGQ = 20))
# 
# Data Descriptives:
#   Number of Observations: 1907
# Number of Groups: 294 
# 
# Model:
#   family: binomial
# link: logit 
# 
# Fit statistics:
#   log.Lik      AIC      BIC
# -623.7227 1257.445 1275.863
# 
# Random effects covariance matrix:
#   StdDev
# (Intercept) 4.022086
# 
# Fixed effects:
#             Estimate Std.Err z-value    p-value
# (Intercept)  -1.6365  0.4293 -3.8119 0.00013792
# time         -0.4045  0.0460 -8.7969    < 1e-04
# treatn       -0.1323  0.5794 -0.2284 0.81934260
# trt.time     -0.1586  0.0719 -2.2055 0.02741885
# 
# Integration:
#   method: adaptive Gauss-Hermite quadrature rule
# quadrature points: 20
# 
# Optimization:
#   method: hybrid EM and quasi-Newton
# converged: TRUE 
#              (Intercept)         time       treatn      trt.time          D_11
# (Intercept)  0.184315166 -0.003454991 -0.162264085  0.0061987920 -0.0128675924
# time        -0.003454991  0.002114794  0.005691257 -0.0018365369 -0.0013049413
# treatn      -0.162264085  0.005691257  0.335745140 -0.0125202414 -0.0023330294
# trt.time     0.006198792 -0.001836537 -0.012520241  0.0051726715 -0.0005864624
# D_11        -0.012867592 -0.001304941 -0.002333029 -0.0005864624  0.0088701084


#VOTING TURNOUT DATA

library(pscl)
#Green and Vavreck (2008) implemented a cluster-randomized experimental design in assessing the
#effects of a voter mobilization treatment in the 2004 U.S. Presidential election. The clusters in this
#design are geographic areas served by a single cable television system.
data("RockTheVote")
dat2 = RockTheVote
print(dat2[1:10,])
res2 = mixed_model(cbind(r,n-r) ~ treated, ~ 1|strata, data=dat2, family=binomial(),
  control=list(nAGQ=20))
print(summary(res2))
Omega2 = solve(res2$Hessian)
print(Omega2)

# Call:
# mixed_model(fixed = cbind(r, n - r) ~ treated, random = ~1 | 
#     strata, data = dat2, family = binomial(), control = list(nAGQ = 20))
# 
# Data Descriptives:
# Number of Observations: 85
# Number of Groups: 40 
# 
# Model:
#  family: binomial
#  link: logit 
# 
# Fit statistics:
#    log.Lik      AIC      BIC
#  -432.0903 870.1807 875.2473
# 
# Random effects covariance matrix:
#                StdDev
# (Intercept) 0.3567023
# 
# Fixed effects:
#             Estimate Std.Err z-value    p-value
# (Intercept)   0.0831  0.0600  1.3865 0.16560042
# treated       0.0990  0.0291  3.4060 0.00065931
# 
# Integration:
# method: adaptive Gauss-Hermite quadrature rule
# quadrature points: 20
# 
# Optimization:
# method: hybrid EM and quasi-Newton
# converged: TRUE 
#               (Intercept)       treated          D_11
# (Intercept)  3.596176e-03 -3.967698e-04 -9.258467e-06
# treated     -3.967698e-04  8.451787e-04  3.926042e-05
# D_11        -9.258467e-06  3.926042e-05  1.441689e-02

bh = res2$post_modes
post_Vars = as.vector(res2$post_vars)
pvar = rep(0,40)
for (i in 1:40) {
  pvar[i] = post_Vars[[i]]
}
cihw = qnorm(0.975)*sqrt(pvar)
lower = bh - cihw
upper = bh + cihw
plot(1, type = "n", xlab = "id", ylab = "bh",
     xlim = c(0, 41), ylim = c(-1.2, 1.2))
id = 1:40
for (i in 1:40) {
  points(id[i],bh[i])
  lines(c(id[i],id[i]),c(lower[i],upper[i]))
}

#TOENAIL EXAMPLE REVISITED
res3a = mixed_model(y ~ time + treatn + trt.time, ~ 1 + time + trt.time |idnew, 
  data=dat1, family=binomial(), control=list(nAGQ=1))
res3b = mixed_model(y ~ time + treatn + trt.time, ~ 1 + time + trt.time |idnew, 
  data=dat1, family=binomial(),
  initial_values = list(betas=res3a$coefficients, D=res3a$D),
  control=list(nAGQ=20))
print(summary(res3b))
Omega3 = solve(res3b$Hessian)
print(round(Omega3,digits=4))
res3c = mixed_model(y ~ time + treatn + trt.time, ~ 1 + time |idnew, 
  data=dat1, family=binomial(), control=list(nAGQ=20))
print(summary(res3c))
Omega3new = solve(res3c$Hessian)
print(round(Omega3new,digits=4))

# > print(summary(res3b))
# 
# Call:
#   mixed_model(fixed = y ~ time + treatn + trt.time, random = ~1 + 
#                 time + trt.time | idnew, data = dat1, family = binomial(), 
#               initial_values = list(betas = res3a$coefficients, D = res3a$D), 
#               control = list(nAGQ = 20))
# 
# Data Descriptives:
#   Number of Observations: 1907
# Number of Groups: 294 
# 
# Model:
#   family: binomial
# link: logit 
# 
# Fit statistics:
#   log.Lik      AIC     BIC
# -547.2523 1114.505 1151.34
# 
# Random effects covariance matrix:
#             StdDev    Corr        
# (Intercept)  8.7304  (Intr)    time
# time         0.9326 -0.6252        
# trt.time     0.3317 -0.3903  0.8438
# 
# Fixed effects:
#             Estimate Std.Err z-value   p-value
# (Intercept)  -2.8766  0.9994 -2.8784 0.0039971
# time         -0.7473  0.2593 -2.8824 0.0039469
# treatn        0.3307  1.2959  0.2552 0.7985820
# trt.time     -0.6590  0.4758 -1.3852 0.1659847
# 
# Integration:
#   method: adaptive Gauss-Hermite quadrature rule
# quadrature points: 20

# Optimization:
#   method: hybrid EM and quasi-Newton
# converged: TRUE 
# > 
#   
#   
# > print(summary(res3c))
# 
# Call:
#   mixed_model(fixed = y ~ time + treatn + trt.time, random = ~1 + 
#                 time | idnew, data = dat1, family = binomial(), control = list(nAGQ = 20))
# 
# Data Descriptives:
#   Number of Observations: 1907
# Number of Groups: 294 
# 
# Model:
#   family: binomial
# link: logit 
# 
# Fit statistics:
#   log.Lik      AIC      BIC
# -548.3806 1110.761 1136.546
# 
# Random effects covariance matrix:
#   StdDev    Corr
# (Intercept)  8.6176        
# time         1.0631 -0.5863
# 
# Fixed effects:
#             Estimate Std.Err z-value  p-value
# (Intercept)  -2.5748  1.2336 -2.0873 0.036865
# time         -0.8798  0.2149 -4.0946  < 1e-04
# treatn       -0.2389  1.3857 -0.1724 0.863136
# trt.time     -0.3665  0.2094 -1.7499 0.080131
# 
# Integration:
#   method: adaptive Gauss-Hermite quadrature rule
# quadrature points: 20
# 
# Optimization:
#   method: hybrid EM and quasi-Newton
# converged: TRUE   

# > Omega3new = solve(res3c$Hessian)
# > print(round(Omega3new,digits=4))
#            (Intercept)    time  treatn trt.time    D_11    D_12    D_22
# (Intercept)      1.5218 -0.0445 -1.1497   0.0824 -0.1031  0.0666 -0.0595
# time            -0.0445  0.0462  0.0631  -0.0222 -0.0102  0.0015 -0.0257
# treatn          -1.1497  0.0631  1.9203  -0.1658  0.0211 -0.0108  0.0243
# trt.time         0.0824 -0.0222 -0.1658   0.0439  0.0004 -0.0006  0.0009
# D_11            -0.1031 -0.0102  0.0211   0.0004  0.0234 -0.0168  0.0152
# D_12             0.0666  0.0015 -0.0108  -0.0006 -0.0168  0.0242 -0.0078
# D_22            -0.0595 -0.0257  0.0243   0.0009  0.0152 -0.0078  0.0284
# > 
