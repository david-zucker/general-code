#STATISTICS FOR DATA SCIENCE - TARGIL 6 - QUESTION 3

### PART (a) ##################################################################

set.seed(3001172)

n = 300
x = rgamma(n,3,0.5)

xgrid = seq(0,15,0.02)
ngrid = length(xgrid)
dentru = dgamma(xgrid,3,0.5)

adj = n^(-0.05)

dd = density(x,kernel='epa',bw='SJ',adjust=adj,from=0,to=15,n=ngrid)
h = dd$bw
plot(dd$x,dd$y,type='l')

k22 = 0.6
varval = k22*dd$y / (n*h)
se = sqrt(varval)
zcrit = qnorm(0.975)
cihw = zcrit*se
f_lo = dd$y - cihw
f_hi = dd$y + cihw
plot(dd$x,dd$y,type='l',xlab='x',ylab='f',ylim=c(0,0.20))
lines(dd$x,f_lo,col='red')
lines(dd$x,f_hi,col='red')
lines(xgrid,dentru,col='green')

### PART (b) ##################################################################

densband = function(x,adjval,nboot) {
  
  n = length(x)
  
  #DENSITY ESTIMATE ON ORIGINAL DATA
  dd_orig = density(x,kernel='epa',bw='SJ',adjust=adjval,from=0,to=15,n=ngrid)
  
  #BOOTSTRAP PROCEDURE
  maxvec = rep(0,nboot)
  for (b in 1:nboot) {
    xboot = sample(x,n,replace=TRUE)
    dd_boot = density(xboot,kernel='epa',bw='SJ',adjust=adjval,from=0,to=15,n=ngrid)
    maxvec[b] = max(abs(dd_boot$y-dd_orig$y))
  }  
  
  #RETURN RESULT
  rslt = quantile(maxvec,probs=0.95,type=4)
  return(rslt)
  
}

dval = densband(x,adj,1000)
f_lo_1 = dd$y - dval
f_hi_1 = dd$y + dval
plot(dd$x,dd$y,type='l',xlab='x',ylab='f',ylim=c(0,0.20))
lines(dd$x,f_lo_1,col='red')
lines(dd$x,f_hi_1,col='red')
lines(xgrid,dentru,col='green')

### PART (c) ##################################################################

dtest = function(settings) {
  n = settings[1]
  adjdel = settings[2]
  nrep = settings[3]
  nboot = settings[4]
  adjval = n^(-adjdel)
  cover = 0
  cihw = rep(0,nrep)
  for (irep in 1:nrep) {
    if (irep %% 10 == 0) print(irep)
    xcur = rgamma(n,3,0.5)
    ddcur = density(xcur,kernel='epa',bw='SJ',adjust=adjval,from=0,to=15,n=ngrid)
    dval = densband(xcur,adjval,nboot)
    covcur = (max(abs(ddcur$y-dentru)) <= dval)
    cover = cover + covcur
    cihw[irep] = dval
  }
  covpr = cover/nrep
  mean_cihw = round(mean(cihw),digits=4)
  rslt = c(n,adjdel,nrep,nboot,covpr,mean_cihw)
  return(rslt)
}

rsltdf = NULL
dt = dtest(c(150,0.025,1000,1000))
rsltdf = rbind(rsltdf,dt)
dt = dtest(c(150,0.05,1000,1000))
rsltdf = rbind(rsltdf,dt)
dt = dtest(c(300,0.025,1000,1000))
rsltdf = rbind(rsltdf,dt)
dt = dtest(c(300,0.05,1000,1000))
rsltdf = rbind(rsltdf,dt)
rsltdf = as.data.frame(rsltdf)
rownames(rsltdf) = NULL
colnames(rsltdf) = c('n', 'adjdel', 'nrep', 'nboot', 'covpr', 'mean_cihw')
print(rsltdf)
  
#     n adjdel nrep nboot covpr mean_cihw
# 1 150  0.025 1000  1000 0.997    0.0534
# 2 150  0.050 1000  1000 0.998    0.0582
# 3 300  0.025 1000  1000 0.991    0.0373
# 4 300  0.050 1000  1000 0.997    0.0411
  
  
  
  



  
  
  



