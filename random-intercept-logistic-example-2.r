#VOTING TURNOUT DATA

library(GLMMadaptive)

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
