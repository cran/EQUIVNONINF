postmys <- function(n,dq,sd,eps1,eps2,tol)   {

ny <- n-1
sigupp <- sd*(ny/qchisq(tol,ny))**0.5
a <- sigupp/2
B <- sigupp

C <- 0.5*C*B;
G <- 0.5*G;

y <- 0

for (j in 1:48)
   { 
     cj <- C[j]
     sigr <- a+cj
     normr2 <- pnorm(sqrt(n)*(eps2*sigr - dq)/sigr)
     normr1 <- pnorm(sqrt(n)*(-eps1*sigr - dq)/sigr)
     normr  <- normr2-normr1
     if (normr != 0)
       { frlog <- log(normr) - n*log(sigr) - 0.5*ny*(sd/sigr)**2 - lgamma(ny/2) +
                  (ny/2)*(log(ny)+2*log(sd)) - (ny/2 - 1)*log(2)
         fr <- exp(frlog)   }
     sigl <- a-cj
     norml2 <- pnorm(sqrt(n)*(eps2*sigl - dq)/sigl)
     norml1 <- pnorm(sqrt(n)*(-eps1*sigl - dq)/sigl)
     norml  <- norml2-norml1
     if (norml != 0)
       { fllog <- log(norml) - n*log(sigl) - 0.5*ny*(sd/sigl)**2 - lgamma(ny/2) +
                  (ny/2)*(log(ny)+2*log(sd)) - (ny/2 - 1)*log(2)
         fl <- exp(fllog)
         y <- y + G[j]*(fr+fl)   }      }

ppost <- B*y

cat(" n =",n," dq =",dq,"   sd =",sd,"   eps1 =",eps1,"   eps2 =",eps2,
    "   tol =",tol,"   PPOST =",ppost)
}
