#STATISTICS FOR DATA SCIENCE - EXERCISE SET 12, QUESTION 3

#SETTINGS
mu = 0.25
sig = 1
n = 20
S = 4

#LOGIT FUNCTION
logit = function(theta) {
  lval = log(theta/(1-theta))
  return(lval)
}

#PRIOR
prior_dist = function(theta,mu,sig) {
  endpt = ((theta==0) | (theta==1))
  rslt = ifelse(endpt,0,
    dnorm(logit(theta),mu,sig)/(theta*(1-theta)))
  return(rslt)
}

#PRIOR TIMES LIKELIHOOD
f1 = function(theta) {
  prior_val = prior_dist(theta,mu,sig)
  lik = (theta^S) * ((1-theta)^(n-S))
  fval = prior_val * lik
  return(fval)
}  

#PRIOR TIMES LIKELIHOOD TIMES theta
#USED FOR COMPUTING POSTERIOR EXPECTATION
f2 = function(theta) {
  prior_val = prior_dist(theta,mu,sig)
  lik = (theta^S) * ((1-theta)^(n-S))
  fval = theta * prior_val * lik
  return(fval)
}  

#NORMALIZING CONSTANT = marginal density of Y evaluated at the data
const_comp = integrate(f1,lower=0, upper=1)
const = const_comp$value

#PART (a)

#POSTERIOR DENSITY
pdens = function(theta) {
  rslt = f1(theta)/const
  return(rslt)
}

thvec = seq(0,1,0.01)
postr = pdens(thvec)
plot(thvec,postr,type='l')

#PREPARE FOR PART (c) - IDENTIFY POSTERIOR MODE
ix = which.max(postr)
thqini = thvec[ix]

#PART (b)

#POSTERIOR CDF
pcdf = function(thval) {
  intcmp = integrate(f1,lower=0,upper=thval)
  prb = intcmp$value/const
  return(prb)
}

prb = pcdf(0.35)
cat('\n')
print(prb)
cat('\n')

#PART (c)

#QUANTILE FUNCTION
qnt = function(q) {
  
  #SET UP FOR NEWTON-RAPHSON ITERATIONS
  maxiter = 50
  tol = 1e-6
  diff = 99
  iter = 0
  thcur = thqini

  #CARRY OUT NEWTON-RAPHSON ITERATIONS
  while ((iter < maxiter) & (diff > tol)) {
    fval = pcdf(thcur) - q
    fdval = pdens(thcur)
    iter = iter + 1
    step = fval/fdval
    thcur = thcur - step
    diff = max(abs(step))
  }
  
  #RETURN ANSWER
  return(thcur)

}

th_lo = qnt(0.05)
th_hi = qnt(0.95)
print(cbind(th_lo,th_hi))
cat('\n')

#PART (d)
nmcmp = integrate(f2,lower=0,upper=1)
thhat = nmcmp$value/const
print(thhat)
cat('\n')

