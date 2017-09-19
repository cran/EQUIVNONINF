mwtie_fr <- function(k,alpha,m,n,eps1_,eps2_,x,y)         {
kpl1 <- k + 1
eqctr <- 0.5 + (eps2_-eps1_)/2
eqleng <- eps1_+eps2_

w <- rep(NA,k)
nx <- rep(NA,k)
ny <- rep(NA,k)
nxc <- rep(NA,kpl1)
nyc <- rep(NA,kpl1)

for (kk in 1:k)
   { w[kk] <- kk
     nx[kk] <- 0
     ny[kk] <- 0
   }


for (i in 1:m)
    for (kk in 1:k)
        if (x[i] == w[kk]) nx[kk] <- nx[kk]+1

for (j in 1:n)
    for (kk in 1:k)
        if (y[j] == w[kk]) ny[kk] <- ny[kk]+1

nxc[1] <- 0
nyc[1] <- 0

for (kk in 2:(k+1)) 
   { nxc[kk] <- nxc[kk-1] + nx[kk-1]
     nyc[kk] <- nyc[kk-1] + ny[kk-1]
   }

wxy <- 0
pih0 <- 0
pihxyy0 <- 0
pihxyy1 <- 0
pihxyy2 <- 0
pihxxy0 <- 0
pihxxy1 <- 0
pihxxy2 <- 0

for (kk in 1:k)
   { pih0 <- pih0 + nx[kk]*ny[kk]
     wxy <- wxy + nx[kk]*nyc[kk]
     pihxyy0 <- pihxyy0 + nx[kk]*ny[kk]*(ny[kk]-1)
     pihxyy1 <- pihxyy1 + nx[kk]*nyc[kk]*(nyc[kk]-1)
     pihxyy2 <- pihxyy2 + nx[kk]*ny[kk]*nyc[kk]
     pihxxy0 <- pihxxy0 + nx[kk]*(nx[kk]-1)*ny[kk]
     pihxxy1 <- pihxxy1 + ny[kk]*(m-nxc[kk+1])*(m-nxc[kk+1]-1)
     pihxxy2 <- pihxxy2 + ny[kk]*nx[kk]*(m-nxc[kk+1])
   }

wxy <- wxy/(m*n)
pih0 <- pih0/(m*n)
wxy_tie <- wxy/(1-pih0)
pihxxy0 <- pihxxy0/(m*(m-1)*n)
pihxyy0 <- pihxyy0/(n*(n-1)*m)
pihxxy1 <- pihxxy1/(m*(m-1)*n)
pihxyy1 <- pihxyy1/(n*(n-1)*m)
pihxxy2 <- pihxxy2/(m*(m-1)*n)
pihxyy2 <- pihxyy2/(n*(n-1)*m)

varhpih0 <- (pih0 - (m+n-1)*pih0**2 + (m-1)*pihxxy0 + (n-1)*pihxyy0) / (m*n)
vargamh <- (wxy - (m+n-1)*wxy**2 + (m-1)*pihxxy1 + (n-1)*pihxyy1) / (m*n)
covh <- ((m-1)*pihxxy2 + (n-1)*pihxyy2 - (m+n-1)*wxy*pih0) / (m*n)
sigmah <- sqrt((1-pih0)**(-2)*vargamh + wxy**2*(1-pih0)**(-4)*varhpih0 + 2*wxy*(1-pih0)**(-3)*covh)
crit <- sqrt(qchisq(alpha,1,(eqleng/2/sigmah)**2))

if (abs((wxy_tie - eqctr)/sigmah) >= crit) rej <- 0 
if (abs((wxy_tie - eqctr)/sigmah) < crit) rej <- 1
if (is.na(sigmah) || is.na(crit)) rej <- 0

cat(" k =",k,"  alpha =",alpha,"  m =",m,"  n =",n,"  eps1_ =",eps1_,"  eps2_ =",eps2_,
    "  WXY_TIE =",wxy_tie,"  SIGMAH =",sigmah,"  CRIT =",crit,"  REJ =",rej)
}


