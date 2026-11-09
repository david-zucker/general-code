#PROGRAM TO DEMONSTRATE SIMULATED STANDARD DEVIATION
#PARAMETER OF INTEREST IS LOG VARIANCE

#gamma distribution with Wasserman parameterization

#PRELIMINARIES
set.seed(342345)
n = 100
alf = 3
bet = 0.5
R = 1000

#ASYMPTOTIC THEORY
M = rep(0,4)
M[1] = alf*bet
for (k in 2:4) {
  M[k] = M[k-1]*(alf+k-1)*bet
}
th = M[2] - M[1]^2
var_th = M[4] - 4*M[1]*M[3] + 8*(M[1]^2)*M[2] - M[2]^2 - 4*(M[1]^4)
var_th = var_th/n
sd_th_asy = sqrt(var_th)
sd_log_th_asy = sd_th_asy/th 

#SIMULATION
lthhat = rep(0,R)
for (r in 1:R) {
  ycur = rgamma(n,alf,scale=bet)
  lthh = log(var(ycur))
  lthhat[r] = lthh
 }
sd_log_th_sim = sd(lthhat)

print(sd_th_asy)
print(sd_log_th_sim)

