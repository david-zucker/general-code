# program to illustrate the drawing of beta(alf,bet) variates
# using rejection sampling

# setup
alf = 2
bet = 3
c = (alf-1)^(alf-1) * (bet-1)^(bet-1) / (alf+bet-2)^(alf+bet-2)
set.seed(24751)

# plot the beta density
x = c(0:100)/100
fx = dbeta(x,alf,bet)
plot(x,fx,type='l')

# draw the sample
nobs = 1000
z = NULL
for (i in 1:nobs) {
  flg = 0
  while (flg == 0) {
    w = runif(1)
    u = runif(1)
    bound = w^(alf-1) * (1-w)^(bet-1) / c
    if (u <= bound) {
      z = c(z,w)
      flg = 1
     }
   }
 }

# estimate the density based on the sample and plot
dens=density(z)
points(dens$x,dens$y)


   
