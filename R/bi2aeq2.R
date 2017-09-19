bi2aeq2 <- function(rho1,rho2,alpha,p1,p2,beta,qlambd)       {

eps <- rho2-1
rho1l <- log(rho1)
rho2l <- log(rho2)
p1l <- log(p1)
q1l <- log(1-p1)
p2l <- log(p2)
q2l <- log(1-p2)
rhoal <- p1l+q2l-q1l-p2l

fakl <- rep(0,2000)
ic1 <- rep(0,2000)
ic2 <- rep(0,2000)
gam1 <- rep(0,2000)
gam2 <- rep(0,2000)
powc <- rep(0,2000)
hypa <- rep(0,2000)
hyp01 <- rep(0,2000)
hyp02 <- rep(0,2000)

raus2 <- 0
problem <- 0
fertig <- 0

fakl[1+0] <- 0
for (i in 1:2000)
   { qi <- i
     fakl[1+i] <- fakl[1+i-1] + log(qi)  }

n <- 5

repeat                                                
  {
GT400 <- 0
GT46  <- 0
GT49  <- 0

qm <- qlambd * n
m <- trunc(qm)
if ((m-qm) < 0) m <- m + 1    
ise <- floor(m*p1 + n*p2)
ixl <- max(0,ise-n)
ixu <- min(ise,m)
ixeins <- ceiling(ise*m/(m+n))
hyp01[2+ixu] <- 0
hyp02[2+ixu] <- 0
hypa[2+ixu] <- 0
ooml <- fakl[1+m+n] - fakl[1+ise] - fakl[1+m+n-ise]

for (j in 1:(ixu-ixl+1))
{ ix <- ixu - j
  hl <- fakl[1+m]-fakl[1+ix+1]-fakl[1+m-ix-1]+fakl[1+n]-
        fakl[1+ise-ix-1]-fakl[1+n-ise+ix+1]
  hyp01[2+ix] <- exp(hl + (ix+1)*rho1l - ooml)
  hyp02[2+ix] <- exp(hl + (ix+1)*rho2l - ooml)
  hypa[2+ix] <- exp(hl + (ix+1)*rhoal - ooml)
  hyp01[2+ix] <- hyp01[2+ix] + hyp01[2+ix+1]
  hyp02[2+ix] <- hyp02[2+ix] + hyp02[2+ix+1]
  hypa[2+ix] <- hypa[2+ix] + hypa[2+ix+1]   }

for (ix in ixl:(ixu-1))
{ hyp01[2+ix] <- hyp01[2+ix] / hyp01[2+ixl-1]
  hyp02[2+ix] <- hyp02[2+ix] / hyp02[2+ixl-1]
  hypa[2+ix] <- hypa[2+ix] / hypa[2+ixl-1]  }

hyp01[2+ixl-1] <- 1
hyp02[2+ixl-1] <- 1
hypa[2+ixl-1] <- 1

k <- ixeins

if (2*k < ise || rho1 != 1/rho2)
  GT400 <- 1
else
{ hrho1k <- hyp01[2+k] - hyp01[2+k+1]
  if (hrho1k >= alpha) GT46 <- 1       }

if (GT46 == 0)
{
  k1 <- min(k+5,ise,m)                        
  repeat                                      
  {
    GT45 <- 0
    k2 <- k1-1
    hrho1c1 <- hyp01[2+k1-1]
    hrho2c1 <- hyp02[2+k1-1]
    
    repeat                                    
    {
      alpha1 <- hrho1c1 - hyp01[2+k2]
      alpha2 <- hrho2c1 - hyp02[2+k2]
      
      if (max(alpha1,alpha2) - alpha > 0) break   
      k2 <- k2 + 1
      if (k2 > ise)
      { GT45 <- 1
        break      }               }          
    
    if (GT45 == 1)
      { k1 <- k1-1
        next        }              
    
                                              
      { k2 <- k2-1
        if (k2 < k1)
          { inc_l <- 1                        
            inc_r <- 1 }
        else
          { k1 <- k1 + 1
            inc_l <- 0
            inc_r <- 1  }   }

      repeat                                  
       {
         bed1 <- 0
         bed2 <- 0
         bed3 <- 0
         alpha1 <- hyp01[2+k1-1] - hyp01[2+k2]
         alpha2 <- hyp02[2+k1-1] - hyp02[2+k2]
         delalph1 <- alpha - alpha1
         delalph2 <- alpha - alpha2
         exhyp11 <- hyp01[2+k1-2] - hyp01[2+k1-1]
         exhyp12 <- hyp01[2+k2] - hyp01[2+k2+1]
         exhyp21 <- hyp02[2+k1-2] - hyp02[2+k1-1]
         exhyp22 <- hyp02[2+k2] - hyp02[2+k2+1]
         dete <- exhyp11*exhyp22 - exhyp12*exhyp21
         gamma1 <- (exhyp22*delalph1-exhyp12*delalph2)/dete
         gamma2 <- (exhyp11*delalph2-exhyp21*delalph1)/dete
         
         if ( (min(gamma1,gamma2) < 0 || max(gamma1,gamma2) >= 1) && 
              inc_l == 0 && inc_r == 1 )
         { k1 <- k1 - 1
           k2 <- k2 - 1
           inc_l <- 1
           inc_r <- 0
           bed1 <- 1     }
         
         if ( (min(gamma1,gamma2) < 0 || max(gamma1,gamma2) >= 1) && 
              inc_l == 1 && inc_r == 0 && bed1 == 0)
         { k2 <- k2 + 1
           inc_l <- 1
           inc_r <- 1
           bed2 <- 1     }
         
         if ( (min(gamma1,gamma2) < 0 || max(gamma1,gamma2) >= 1) && 
              inc_l == 1 && inc_r == 1 && bed2 == 0)
         { k1 <- k1 - 1
           bed3 <- 1
           break         }                     
         
         if ( bed1 == 0 && bed2 == 0 && bed3 == 0)
         { GT49 <- 1
         break       }              }           
         
       if (GT49 == 1) break                     
    
      
          } 
                              }         

if (GT46 == 1)                                      
{ icl <- k
  icr <- k
  gamma1 <- alpha/hrho1k
  gamma2 <- gamma1         }

ic1e <- k1-1                             
ic2e <- k2+1
pagtc1 <- hypa[2+ic1e]
pagtc1mi <- hypa[2+ic1e-1]
pagtc2mi <- hypa[2+ic2e-1]
pagtc2 <- hypa[2+ic2e]
ic1eqc2 <- 0
if (ic1e == ic2e) ic1eqc2 <- 1
pownre <- (pagtc1 - pagtc2mi) * (1-ic1eqc2)
powe <- pownre + gamma1 * (pagtc1mi-pagtc1) * (1-0.5*ic1eqc2) +
        gamma2 * (pagtc2mi-pagtc2) * (1-0.5*ic1eqc2)

if (powe-beta > 0)
  { raus2 <- 1
    break }                                   

if (powe-beta <= 0)
  { n <- n + 5
    if ((1+qlambd)*n-2000 >= 0)
      { cat("Summe der benoetigten Stichprobenumfaenge > 2000")
        break   }  }                               

                                  }      

if (raus2 == 1)
{ pow <- pwexacta(m,n,alpha,q1l,q2l,p1l,p2l,rho1,rho2,rhoal,fakl)
  
  if (pow - beta < 0)
  {
    repeat   
    {
      n <- n + 5                          
      if ((1+qlambd)*n-2000 >= 0)
      { cat("Summe der benoetigten Stichprobenumfaenge > 2000")
        problem <- 1
        break   }
      
      qm <- qlambd*n                       
      m <- trunc(qm)
      if ((m-qm) < 0) m <- m + 1
      pow <- pwexacta(m,n,alpha,q1l,q2l,p1l,p2l,rho1,rho2,rhoal,fakl)
      if (pow - beta >= 0) break
    }
    if (problem == 0)
    {
      repeat
      {
        n <- n - 1                                
        qm <- qlambd*n
        m <- trunc(qm)
        if ((m-qm) < 0) m <- m + 1
        pow <- pwexacta(m,n,alpha,q1l,q2l,p1l,p2l,rho1,rho2,rhoal,fakl)
        if (pow - beta < 0) break      }
    
    n <- n + 1                                
    qm <- qlambd*n
    m <- trunc(qm)
    if ((m-qm) < 0) m <- m + 1
    pow <- pwexacta(m,n,alpha,q1l,q2l,p1l,p2l,rho1,rho2,rhoal,fakl)
    fertig <- 1                     }  }
  
  if (pow - beta == 0) fertig <- 1
  
  if (pow - beta > 0) 
  {
    repeat
    {
      n <- n - 5                                 
      qm <- qlambd*n
      m <- trunc(qm)
      if ((m-qm) < 0) m <- m + 1
      pow <- pwexacta(m,n,alpha,q1l,q2l,p1l,p2l,rho1,rho2,rhoal,fakl)
      if (pow - beta < 0) break      }
    
    repeat
    {
      n <- n + 1                                
      qm <- qlambd*n
      m <- trunc(qm)
      if ((m-qm) < 0) m <- m + 1
      pow <- pwexacta(m,n,alpha,q1l,q2l,p1l,p2l,rho1,rho2,rhoal,fakl)
      fertig <- 1 
      if (pow - beta >= 0) break  }        }
  
  cat("\n"," EPS =",eps,"   ALPHA =",alpha,"\n",
      " P1 =",p1,"   P2 =",p2,"\n",
      " BETA =",beta,"   LAMBDA =",qlambd,"\n","\n",
      " M = ",m,"  N = ",n,"     POW = ",pow)        
  
  }                                   
                                                          }
  



