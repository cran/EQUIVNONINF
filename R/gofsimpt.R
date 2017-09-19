gofsimpt <- function(alpha,n,k,eps,x,pio)                {
pih <- rep(NA,k)

for (J in 1:k) pih[J] <- x[J] / n

dsqpih_0 <- 0
vnsq_1 <- 0

for (J in 1:k)
   { dsqpih_0 <- dsqpih_0 + (pih[J] - pio[J])**2
     vnsq_1 <- vnsq_1 + (pih[J] - pio[J])**2 * pih[J]
   }

vnsq_2 <- 0

for (J1 in 1:k)
    for (J2 in 1:k)
        vnsq_2 <- vnsq_2 + (pih[J1]-pio[J1])*(pih[J2]-pio[J2])*pih[J1]*pih[J2]

vnsq <- (4/n) * (vnsq_1 - vnsq_2)
vn_n <- sqrt(vnsq)
epsaksq <- eps**2
crit <- epsaksq - qnorm(1-alpha)*vn_n
rej <- 0

if (is.na(dsqpih_0) == FALSE && dsqpih_0 < crit)  rej <- 1

cat("  n =",n," alpha =",alpha,"  eps =",eps,"   x(1,K) =",x,"  pio(1,K) =",pio,
    "   DSQPIH_0 =",dsqpih_0,"  VN_N =",vn_n,"  CRIT =",crit,"  REJ =",rej)
}


 

