zcrit = qnorm(0.975)
expit = function(u) exp(u)/(1+exp(u))
a = c(1, 50, 0, 200, 100, 1)
lr1 = glm(died10 ~ age + diab + chol + map + smoke, family=binomial(link=logit))
print(summary(lr1))
lammat = vcov(lr1)
print(lammat)
beta = coef(lr1)
psihat = sum(a*beta)
print(psihat)
phat = expit(psihat)
se_psi = sqrt(t(a) %*% lammat %*% a)
cihw_psi = zcrit*se_psi
psi_lo = psihat - cihw_psi
psi_hi = psihat + cihw_psi
p_lo = expit(psi_lo)
p_hi = expit(psi_hi)
print(cbind(p_lo,p_hi))