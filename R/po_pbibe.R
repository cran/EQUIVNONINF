po_pbibe <- function(n,eps,pio,zq,s,tol,sw,ihmax)   {

epstil <- log(1+eps)
sigeps0 <- epstil/qnorm((1+pio)/2)
a <- sigeps0/2
B <- sigeps0

sigl <- rep(NA,48)
sigr <- rep(NA,48)
ql <- rep(NA,48)
qr <- rep(NA,48)

      C <- 0.5*C*B;
      G <- 0.5*G;


for (j in 1:48)
   {
   sigl[j] <-  a - C[j]
   sigr[j] <-  a + C[j]
   for (ny in 1:2)
     {
      if (ny == 1) epst_sig <- epstil/sigl[j]
      if (ny == 2) epst_sig <- epstil/sigr[j]
      psi <- 0
      targ <- 2*pnorm(epst_sig) - 1
      while (targ >= pio)
           { psi <- psi + sw
             targ <- pnorm(epst_sig-psi) - pnorm(-epst_sig-psi)   }
      psi1 <- psi-sw
      psi2 <- psi
      ih <- 0
      repeat
        {
          psi0 <- (psi1+psi2)/2
          ih <- ih + 1
          targ <- pnorm(epst_sig-psi0) - pnorm(-epst_sig-psi0)
          if (abs(targ-pio) < tol || ih >= ihmax) break  else
          { if (targ <= pio-tol) psi2 <- psi0
            if (targ >= pio+tol) psi1 <- psi0  }
        }
      if (ny == 1) ql[j] <- psi0
      if (ny == 2) qr[j] <- psi0      }  }

lghn_1 <- lgamma((n-1)/2)
y <- 0

for (j in 1:48)
   { phidffl <- pnorm(sqrt(n)*( ql[j]-zq/sigl[j])) -
                pnorm(sqrt(n)*(-ql[j]-zq/sigl[j]))
     vaul <- sqrt(n-1)*s/sigl[j]
     ldsigl <- (n-1)*log(vaul) - vaul**2/2 - (n-3)/2*log(2) - lghn_1 - log(sigl[j])
     igdsigl <- phidffl*exp(ldsigl)
     
     phidffr <- pnorm(sqrt(n)*( qr[j]-zq/sigr[j])) -
                pnorm(sqrt(n)*(-qr[j]-zq/sigr[j]))
     vaur <- sqrt(n-1)*s/sigr[j]
     ldsigr <- (n-1)*log(vaur) - vaur**2/2 - (n-3)/2*log(2) - lghn_1 - log(sigr[j])
     igdsigr <- phidffr*exp(ldsigr)
     y <- y + G[j]*(igdsigl + igdsigr)    }

po_pbibe <- B*y

cat(" n =",n," eps =",eps,"   pio =",pio,"   zq =",zq,"   s =",s,"   tol =",tol,
    "   sw =",sw,"   IHMAX =",ihmax,"   PO_PBIBE =",po_pbibe)  
}



