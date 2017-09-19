mawi <- function(alpha,m,n,eps1_,eps2_,x,y)                {
eqctr <- 0.5 + (eps2_-eps1_)/2 
eqleng <- eps1_ + eps2_

wxy <- 0
pihxxy <- 0
pihxyy <- 0

for (i in 1:m)
    for (j in 1:n)
        wxy <- wxy + trunc(0.5*(sign(x[i] - y[j]) + 1))

for (i in 1:m)
    for (j1 in 1:(n-1))
        for (j2 in (j1+1):n)
            pihxyy <- pihxyy + trunc(0.5*(sign(x[i] - max(y[j1],y[j2])) + 1))

for (i1 in 1:(m-1))
    for (i2 in (i1+1):m)
        for (j in 1:n)
            pihxxy <- pihxxy + trunc(0.5*(sign(min(x[i1],x[i2]) - y[j]) + 1))

wxy <- wxy / (m*n)
pihxxy <- pihxxy*2 / (m*(m-1)*n)
pihxyy <- pihxyy*2 / (n*(n-1)*m)
sigmah <- sqrt((wxy-(m+n-1)*wxy**2+(m-1)*pihxxy+(n-1)*pihxyy)/(m*n))

crit <- sqrt(qchisq(alpha,1,(eqleng/2/sigmah)**2))

if (abs((wxy-eqctr)/sigmah) >= crit) rej <- 0
if (abs((wxy-eqctr)/sigmah) < crit)  rej <- 1

if (is.na(sigmah) || is.na(crit)) rej <- 0

cat(" alpha =",alpha,"  m =",m,"  n =",n,"  eps1_ =",eps1_,"  eps2_ =",eps2_,
 "\n","W+ =",wxy,"  SIGMAH =",sigmah,"  CRIT =",crit,"  REJ =",rej)
}
