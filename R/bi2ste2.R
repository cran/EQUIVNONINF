bi2ste2 <-
function(eps,alpha,p1,p2,bet,qlambd)       {

raus2 <- 0
problem <- 0
fertig <- 0

rho0l <- log(1-eps)
p1l <- log(p1)
q1l <- log(1-p1)
p2l <- log(p2)
q2l <- log(1-p2)
rhoal <- p1l + q2l - q1l - p2l

fakl <- rep(NA,2000)
hyp0 <- rep(NA,2000)
hypa <- rep(NA,2000)

fakl[1+0] <- 0
for (i in 1:2000)
   { qi <- i
     fakl[1+i] <- fakl[1+i-1] + log(qi) }

n <- 5
repeat                       
{
qm <- qlambd * n
m <- trunc(qm)
if ((m-qm) < 0) m <- m + 1    
ise <- floor(m*p1+n*p2)
ixl <- max(0,ise-n)
ixu <- min(ise,m)
hyp0[2+ixu] <- 0
hypa[2+ixu] <- 0
ooml <- fakl[1+m+n] - fakl[1+ise] - fakl[1+m+n-ise]

for (j in 1:(ixu-ixl+1))
   { ix <- ixu - j
     hl <- fakl[1+m]-fakl[1+ix+1]-fakl[1+m-ix-1]+fakl[1+n]-
           fakl[1+ise-ix-1]-fakl[1+n-ise+ix+1]
     hyp0[2+ix] <- exp(hl + (ix+1)*rho0l - ooml)
     hypa[2+ix] <- exp(hl + (ix+1)*rhoal - ooml)
     hyp0[2+ix] <- hyp0[2+ix] + hyp0[2+ix+1]
     hypa[2+ix] <- hypa[2+ix] + hypa[2+ix+1]   }
for (ix in ixl:(ixu-1))
   { hyp0[2+ix] <- hyp0[2+ix] / hyp0[2+ixl-1]
     hypa[2+ix] <- hypa[2+ix] / hypa[2+ixl-1]  }

hyp0[2+ixl-1] <- 1
hypa[2+ixl-1] <- 1
k <- ixu+1
k <- k - 1
size <- hyp0[2+k]
while (size-alpha <= 0)
    { k <- k - 1
      size <- hyp0[2+k]  }
ice <- k + 1
game <- (alpha-hyp0[2+ice]) / (hyp0[2+k]-hyp0[2+ice])
pow <- hypa[2+ice]*(1-game) + game*hypa[2+k]

if (pow-bet > 0)
  { raus2 <- 1
    break }          
if (pow-bet <= 0)
    { n <- n + 5
      if ((1+qlambd)*n-2000 >= 0)
        { cat("Summe der benoetigten Stichprobenumfaenge > 2000")
          break   }  }       }         # Ende repeat  

if (raus2 == 1)
{
pow <- pwexact(m,n,eps,alpha,p1,p2)                       

if (pow - bet < 0)
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
     pow <- pwexact(m,n,eps,alpha,p1,p2)
     if (pow - bet >= 0) break
                                      }
if (problem == 0)
  {
  repeat
{
n <- n - 1                                
qm <- qlambd*n
m <- trunc(qm)
if ((m-qm) < 0) m <- m + 1
pow <- pwexact(m,n,eps,alpha,p1,p2)
if (pow - bet < 0) break      }

n <- n + 1                                
qm <- qlambd*n
m <- trunc(qm)
if ((m-qm) < 0) m <- m + 1
pow <- pwexact(m,n,eps,alpha,p1,p2)
fertig <- 1                     }   }

if (pow - bet == 0) fertig <- 1

if (pow - bet > 0) 
{
repeat
{
n <- n - 5                                 
qm <- qlambd*n
m <- trunc(qm)
if ((m-qm) < 0) m <- m + 1
pow <- pwexact(m,n,eps,alpha,p1,p2)
if (pow - bet < 0) break      }

repeat
{
n <- n + 1                                
qm <- qlambd*n
m <- trunc(qm)
if ((m-qm) < 0) m <- m + 1
pow <- pwexact(m,n,eps,alpha,p1,p2)
fertig <- 1 
if (pow - bet >= 0) break  }        }     

cat(" EPS = ",eps,"   ALPHA = ",alpha,"\n",
    "P1 = ",p1,"   P2 = ",p2,"\n",
    "BETA = ",bet,"   LAMBDA = ",qlambd,"\n","\n",
    "M = ",m,"   N = ",n,"        POW = ",pow)

erg <- c(m,n,pow)

                                         }  


}
