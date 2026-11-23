#STATISTICS FOR DATA SCIENCE - EXERCISE SET 5, QUESTION 3

set.seed(1458148)
options(max.print=100000)

#READ DATA
indat = read.table('bmi.txt')
xorig = indat[,1]

#remove the extreme outlier for purposes of this demonstration
ix = which(xorig>40)
xorig = xorig[-ix]

#SETUPS
n = length(xorig)
R = 1000
M = 100

#FUNCTION TO READ DATA AND COMPUTE SKEWNESS
skewcomp = function(xcur) {
  mm = rep(0,3)
  for (j in 1:3) {
    mm[j] = mean(xcur^j)
  }
  var = mm[2] - mm[1]^2
  sig = sqrt(var)
  sig3 = sig^3
  mu3hat = mm[3] - 3*mm[2]*mm[1] + 2*(mm[1]^3)
  gamhat = mu3hat/sig3
  return(gamhat)
} 

#MAIN PROGRAM

#COMPUTE ESTIMATE BASED ON ORIGINAL SAMPLE
gamhat_orig = skewcomp(xorig)

#DOUBLE BOOTSTRAP
gam_hat_star = rep(0,R)
gam_hat_dbl_star = matrix(0,R,M)
#begin outer loop
for (r in 1:R) {
  if ((r %% 50) == 0) print(noquote(paste('r =',r)))
  xboot = sample(xorig, n, replace=TRUE)
  gam_hat_star[r] = skewcomp(xboot)
  #begin inner loop
  for (m in 1:M) {
    xsam_bt_bt = sample(xboot, n, replace=TRUE)
    gam_hat_dbl_star[r,m] = skewcomp(xsam_bt_bt)
  }
  #end inner loop
}
#end outer loop

#FINAL COMPUTATIONS
gam_star_bar = mean(gam_hat_star)
gam_st_st_bar = sum(gam_hat_dbl_star)/(R*M)
gamhat_corr = 3*(gamhat_orig - gam_star_bar) + gam_st_st_bar
rslt = as.data.frame(cbind(gamhat_orig,gamhat_corr))
colnames(rslt) = c('gamhat','gamhat_corr')
print(rslt)

