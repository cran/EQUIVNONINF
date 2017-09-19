bi2rlv2 <- function(rho1,rho2,alpha,p1,p2,beta,qlambd)  {

cat("\n"," rho1 =",rho1," rho2 =",rho2," alpha =",alpha," p1 =",p1,"  p2 =",p2,
     " beta =",beta,"  lambda =",qlambd)

alpha <- 1-alpha
beta <- 1-beta
rho1l <- log(rho1)
rho2l <- log(rho2)
p1l <- log(p1)
q1l <- log(1-p1)
p2l <- log(p2)
q2l <- log(1-p2)
rhoal <- p1l+q2l-q1l-p2l

problem <- 0

fakl <- rep(0,2000)

fakl[1+0] <- 0
for (i in 1:2000)
   { qi <- i
     fakl[1+i] <- fakl[1+i-1] + log(qi)  }

n <- 100

repeat                                           
  { 
    qm <- qlambd * n
    m <- trunc(qm)
    if ((m-qm) < 0) m <- m + 1    
    pow <- pwexacta(m,n,alpha,q1l,q2l,p1l,p2l,rho1,rho2,rhoal,fakl)
    if((pow-beta) < 0) break
    n <- n+100
    if ((1+qlambd)*n-2000 > 0)
      { cat("Summe der benoetigten Stichprobenumfaenge > 2000")
        problem <- 1
        break   }
                                               }               

if (problem == 0)                                      
  
  { n <- n-100
    repeat                                          
      { n <- n+10
        qm <- qlambd * n
        m <- trunc(qm)
        if ((m-qm) < 0) m <- m + 1    
        pow <- pwexacta(m,n,alpha,q1l,q2l,p1l,p2l,rho1,rho2,rhoal,fakl)
        if((pow-beta) < 0) break
      
                                                  }  
  
  n <- n-10                                    
  repeat                                        
  { n <- n+1
    qm <- qlambd * n
    m <- trunc(qm)
    if ((m-qm) < 0) m <- m + 1   
    pow <- pwexacta(m,n,alpha,q1l,q2l,p1l,p2l,rho1,rho2,rhoal,fakl)
    if((pow-beta) < 0) break
  
                                                 }  
  
  pow <- 1-pow                                
  
  cat("\n","  m = ",m,"  n = ",n,"     POW = ",pow)        
  
  
                                               }    
                                                                            }