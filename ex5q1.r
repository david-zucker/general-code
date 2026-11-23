#STATISTICS FOR DATA SCIENCE - EXERCISE SET 5, QUESTION 1

set.seed(5074283)

#exact calculation
exact = ppois(10,10)
print(exact)

#simulation
nreps = 100000
xbar_vec = rep(0,nreps)
for (i in 1:nreps) {
  datcur = rpois(5,2)
  xbar_vec[i] = mean(datcur)
}
simprb = mean((xbar_vec <= 2))
print(simprb)
print(exact-simprb)


    
  