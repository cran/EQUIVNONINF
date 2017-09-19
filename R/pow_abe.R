pow_abe <- function(m,n,alpha,del_0,del,sig) {

mpln <- m + n
lcon <- 0.5*log(mpln-2) + (2-mpln/2)*log(2) - lgamma(mpln/2 - 1)
sqmn <- sqrt(m*n/mpln)
t_al <- qt(1-alpha,mpln-2)
vm <- sqmn*2*(del_0/sig)/t_al

a <- vm/2
B <- vm

C <- 0.5*C*B;
G <- 0.5*G;

y <- 0

for (j in 1:48)
   {
    cj <- C[j]
    vr <- a + cj
    ur <- sqrt(mpln-2)*vr
    normr2 <- pnorm(sqmn*2*(del_0-del)/sig - vr*t_al)
    normr1 <- pnorm(-sqmn*2*(del_0+del)/sig + vr*t_al)
    normr  <- normr2 - normr1
    if (normr != 0)
      { itgdrlog <- log(normr) + lcon - ur**2/2 + (mpln-3)*log(ur)
        itgdr <- exp(itgdrlog)  }
    vl <- a-cj
    ul <- sqrt(mpln-2)*vl
    norml2 <- pnorm(sqmn*2*(del_0-del)/sig - vl*t_al)
    norml1 <- pnorm(-sqmn*2*(del_0+del)/sig + vl*t_al)
    norml  <- norml2 - norml1
    if (norml != 0)
      { itgdllog <- log(norml) + lcon - ul**2/2 + (mpln-3)*log(ul)
        itgdl <- exp(itgdllog) 
        y <-  y + G[j]*(itgdr+itgdl)  }     }

pow_abe <- B*y

cat(" m =",m,"  n =",n," alpha =",alpha,"   del_0 =",del_0,"   del =",del,
    "   sig =",sig,"   POW_ABE =",pow_abe) 
}


