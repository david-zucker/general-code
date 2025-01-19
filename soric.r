#SELECTIVE INFERENCE EXAMPLE 2

set.seed(9876)

n = 10000
sdth = 0.2
za = qnorm(0.95)

theta = sdth*rnorm(n)
eps = rnorm(n)
z = theta + eps
tlo = z-za
thi = z+za

cover = (theta >= tlo)*(theta <= thi)
zro.out = (tlo > 0) | (thi < 0)

rte1 = mean(cover)
rte2 = mean(cover[which(zro.out)])
p.zro.out = mean(zro.out)

print(cbind(n,sdth,rte1,rte2,p.zro.out))

#         n sdth  rte1      rte2 p.zro.out
#[1,] 10000  0.2 0.896 0.1909594    0.1084