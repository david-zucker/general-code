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
