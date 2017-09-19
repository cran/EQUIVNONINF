bi2aeq1 <- function(m,n,rho1,rho2,alpha,p1,p2)  {

fakl <- rep(0,2000)
hypa <- rep(0,2000)
hyp01 <- rep(0,2000)
hyp02 <- rep(0,2000)

rho1l <- log(rho1)
rho2l <- log(rho2)
p1l <- log(p1)
q1l <- log(1-p1)
p2l <- log(p2)
q2l <- log(1-p2)
rhoal <- p1l+q2l-q1l-p2l

fakl[1+0] <- 0
for (i in 1:2000)
{ qi <- i
  fakl[1+i] <- fakl[1+i-1] + log(qi)  }

nn <- m + n
is <- 0
oom01l <- 0
oom02l <- 0
oomal <- 0
pownr <- 0
powc <- alpha
probs <- exp(m*q1l + n*q2l)
pow <- powc*probs
for (is in 1:(nn-1))                               
   {  
      GT300 <- 0
      GT36  <- 0
      GT39  <- 0
      
      ixl <- max(0,is-n)
      ixu <- min(is,m)
      ixeins <- ceiling(is*m/nn)
      hyp01[2+ixu] <- 0
      hyp02[2+ixu] <- 0
      hypa[2+ixu] <- 0
      for (j in 1:(ixu-ixl+1))                     
         { ix <- ixu-j
           hl <- fakl[1+m]-fakl[1+ix+1]-fakl[1+m-ix-1]+fakl[1+n]-fakl[1+is-ix-1]-fakl[1+n-is+ix+1]
           hyp01[2+ix] <- exp(hl + (ix+1)*rho1l - oom01l)
           hyp02[2+ix] <- exp(hl + (ix+1)*rho2l - oom02l)
           hypa[2+ix] <- exp(hl + (ix+1)*rhoal - oomal)
           hyp01[2+ix] <- hyp01[2+ix] + hyp01[2+ix+1]
           hyp02[2+ix] <- hyp02[2+ix] + hyp02[2+ix+1]
           hypa[2+ix] <- hypa[2+ix] + hypa[2+ix+1]  }      
      oom01l <- oom01l + log(hyp01[2+ixl-1])
      oom02l <- oom02l + log(hyp02[2+ixl-1])
      oomal <- oomal + log(hypa[2+ixl-1])
      for (ix in ixl:(ixu-1))                              
         { hyp01[2+ix] <- hyp01[2+ix] / hyp01[2+ixl-1]
           hyp02[2+ix] <- hyp02[2+ix] / hyp02[2+ixl-1]
           hypa[2+ix] <- hypa[2+ix] / hypa[2+ixl-1]  }     
      hyp01[2+ixl-1] <- 1
      hyp02[2+ixl-1] <- 1
      hypa[2+ixl-1] <- 1
      
      k <- ixeins
      
      if (2*k < is || rho1 != 1/rho2)
        GT300 <- 1
      else
      { hrho1k <- hyp01[2+k] - hyp01[2+k+1]
      if (hrho1k >= alpha) GT36 <- 1       }
      
      if (GT36 == 0)
      {
        k1 <- min(k+5,is,m)                        
        repeat                                      
        {
          GT35 <- 0
          k2 <- k1-1
          hrho1c1 <- hyp01[2+k1-1]
          hrho2c1 <- hyp02[2+k1-1]
          
          repeat                                    
          {
            alpha1 <- hrho1c1 - hyp01[2+k2]
            alpha2 <- hrho2c1 - hyp02[2+k2]
            
            if (max(alpha1,alpha2) - alpha > 0) break    
            k2 <- k2 + 1
            if (k2 > is)
            { GT35 <- 1
              break      }               }          
          
          if (GT35 == 1)
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
              if (k1 == ixl)
              {gamma1 <- 0
               gamma2 <- (alpha - (1-hyp01[2+k2]))/(hyp01[2+k2] - hyp01[2+k2+1]) }  
              if (k2 >= min(is,m) )
              {gamma2 <- 0
               gamma1 <- (alpha - hyp01[2+min(is,m,k1-1)])/
                 (hyp01[2+min(is,m,k1-2)] - hyp01[2+min(is,m,k1-1)] ) } 
              if (k1 > ixl && k2 < min(is,m))
              {delalph1 <- alpha - alpha1
               delalph2 <- alpha - alpha2
               exhyp11 <- hyp01[2+k1-2] - hyp01[2+k1-1]
               exhyp12 <- hyp01[2+k2] - hyp01[2+k2+1]
               exhyp21 <- hyp02[2+k1-2] - hyp02[2+k1-1]
               exhyp22 <- hyp02[2+k2] - hyp02[2+k2+1]
               det <- exhyp11*exhyp22 - exhyp12*exhyp21
               
               gamma1 <- (exhyp22*delalph1-exhyp12*delalph2)/det
               gamma2 <- (exhyp11*delalph2-exhyp21*delalph1)/det  }
              
              if ( (min(gamma1,gamma2) < 0 || max(gamma1,gamma2) >= 1) && 
                   inc_l == 0 && inc_r == 1 )
              { k1 <- k1 - 1
                k2 <- k2 - 1
                inc_l <- 1
                inc_r <- 0
                bed1 <- 1     
                next         }
              
              if ( (min(gamma1,gamma2) < 0 || max(gamma1,gamma2) >= 1) && 
                   inc_l == 1 && inc_r == 0)
              { k2 <- k2 + 1
                inc_l <- 1
                inc_r <- 1
                bed2 <- 1     
                next          }
              
              if ( (min(gamma1,gamma2) < 0 || max(gamma1,gamma2) >= 1) && 
                   inc_l == 1 && inc_r == 1)
              { k1 <- k1 - 1
                bed3 <- 1
                break         }                     
              
              if ( bed1 == 0 && bed2 == 0 && bed3 == 0)
              { GT39 <- 1
                break       }                       
              
                                       }          
          
          if (GT39 == 1) break                     
          
          
                                                    } 
        
        
                                            }          
      
      if (GT36 == 1)                                      
      { c1 <- k
        c2 <- k
        gamma1 <- alpha/hrho1k
        gamma2 <- gamma1         }
      
      ic1 <- k1-1                                     
      ic2 <- k2+1
      gam1 <- gamma1
      gam2 <- gamma2
      pagtc1 <- hypa[2+ic1]
      pagtc1mi <- hypa[2+ic1-1] 
      pagtc2mi <- hypa[2+ic2-1] 
      pagtc2 <- hypa[2+ic2]
      ic1eqc2 <- 0
      if (ic1 == ic2) ic1eqc2 <- 1
      pownrc <- (pagtc1 - pagtc2mi) * (1-ic1eqc2)
      powc <- pownrc + gam1*(pagtc1mi-pagtc1) * (1-0.5*ic1eqc2) +
              gam2* (pagtc2mi-pagtc2) * (1-0.5*ic1eqc2)
      probs <- 0
      for (ix in ixl:ixu)
         { bl <- fakl[1+n]-fakl[1+is-ix]-fakl[1+n-is+ix] + 
                (is-ix)*p2l + (n-is+ix)*q2l + fakl[1+m] - fakl[1+ix] - 
                fakl[1+m-ix] + ix*p1l + (m-ix)*q1l
           probs <- probs + exp(bl)   }
     
      pownr <- pownr + pownrc*probs
      pow <- pow + powc*probs        
      }    #  Ende Do 9

powc <- alpha
probs <- exp(m*p1l + n*p2l)
pow <- pow + powc*probs
cat("\n","m= ",m,"  n= ",n,"rho1= ",rho1,"rho2= ",rho2,"\n","  POWNR=",pownr,"  POW=",pow)   }

                     
