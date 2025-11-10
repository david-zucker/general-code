setwd("C:/Users/owner/Dropbox/WORK/Zucker/Stat for DS/Targilim/5786")

#READ DATA
indat = read.table('bmi.txt')
x = indat[,1]
n = length(x)

#COMPUTE REQUIRED SAMPLE MOMENTS
mm = rep(0,6)
for (j in 1:6) {
  mm[j] = mean(x^j)
}

#QUESTION 5

psihat = mm[2] - mm[1]^2
varpsi = mm[4] - 4*mm[1]*mm[3] + 8*(mm[1]^2)*mm[2] - mm[2]^2 - 4*(mm[1]^4)
varpsi = varpsi/n
sdpsi = sqrt(varpsi)
cvpr = 0.95
zcrit = qnorm(1-(1-cvpr)/2)
cihw1 = zcrit*sdpsi
psi_lo = psihat - cihw1
psi_hi = psihat + cihw1
print(cbind(psihat, cihw1, psi_lo, psi_hi))

#QUESTION 6

eta_hat = log(psihat)
cihw2 = zcrit*sdpsi / psihat
eta_lo = eta_hat - cihw2
eta_hi = eta_hat + cihw2
print(cbind(eta_hat, cihw2, eta_lo, eta_hi))
psi_lo_2 = exp(eta_lo)
psi_hi_2 = exp(eta_hi)
print(cbind(psi_lo_2, psi_hi_2))

#QUESTION 7

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
gamhat = mu3hat/sig^3

G = 3*mm[1]*mu3hat/sig5 + (6*(mm[1]^2) - 3*mm[2])/sig3
H = -(3/2)*mu3hat/sig5 - 3*mm[1]/sig3
K = 1/sig3

tau2 = A*(G^2) + 2*B*G*H + 2*D*G*K + C*(H^2) + 2*E*H*K + F*(K^2)

cihw3 = zcrit * (1/sqrt(n)) * sqrt(tau2)
gam_lo = gamhat - cihw3
gam_hi = gamhat + cihw3
print(cbind(gamhat, cihw3, gam_lo, gam_hi))

