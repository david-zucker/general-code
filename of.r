
#Exercise Set 3, Question 4
#Old Faithful geyser data

setwd("C:/Users/owner/Dropbox/WORK/Zucker/Stat for DS/Targilim/5785")

indat = read.table('ofdat.txt')
colnames(indat) = c('ObsNum', 'eruptions', 'waiting')
n = nrow(indat)

#PART (a)
count.70 = sum(indat$waiting <= 70)
Fn.70 = count.70/n
cnf.lvl = 0.90
alf = 1 - cnf.lvl
zcrit = qnorm(1-alf/2)
cihw = zcrit*sqrt(Fn.70*(1-Fn.70)/n)
F_L = Fn.70 - cihw
F_U = Fn.70 + cihw

#PART (b)

eps = sqrt((1/(2*n))*log(2/alf))

my.edf = function(x) {
  x1 = sort(unique(x))
  n1 = length(x1)
  Fn = rep(0,n1)
  for (j in 1:n1) {
    Fn[j] = mean(x <= x1[j])
  }
  ans = list(x1=x1, Fn=Fn)
  return(ans)
}

Fn.rslt = my.edf(indat$waiting)
plot(Fn.rslt$x1, Fn.rslt$Fn, xlab='x',ylab='F(x)', type='l')
lines(Fn.rslt$x1, Fn.rslt$Fn-eps,col='red')
lines(Fn.rslt$x1, Fn.rslt$Fn+eps,col='red')


