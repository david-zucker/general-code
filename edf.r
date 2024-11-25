#Exercise Set 4, Question 6

#I'm using the fact that the sample points are unique

x = c(0.43, 0.85, 0.94, 0.24, 0.64)

xs = sort(x)
n = length(xs)
cdf = (1:n)/n

plot(xs,cdf, xlim=c(0,1.2), ylim=c(0,1), pch=19)
lines(c(0,xs[1]),c(0,0))
points(c(xs[1]),c(0))
for (i in 2:n) {
  lines(c(xs[i-1],xs[i]), c(cdf[i-1],cdf[i-1]))
  points(c(xs[i]),c(cdf[i-1]))
} 
lines(c(xs[n],1),c(1,1))
  
  



