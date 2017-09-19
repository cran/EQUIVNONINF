bi2ste1 <-
function(m,n,eps,alpha,p1,p2)    {


ic <- rep(NA,2000)
gam <- rep(NA,2000)
pownrc <- rep(NA,2000)
powc <- rep(NA,2000)
fakl <- rep(NA,2000)
hyp0 <- rep(NA,2000)
hypa <- rep(NA,2000)

nn <- m + n
rho0l <- log(1-eps)
p1l <- log(p1)
q1l <- log(1-p1)
p2l <- log(p2)
q2l <- log(1-p2)
rhoal <- p1l + q2l - q1l -p2l

fakl[1+0] <- 0
for (i in 1:2000)
{ qi <- i
fakl[1+i] <- fakl[1+i-1] + log(qi) }

is <- 0
oom0l <- 0
oomal <- 0
ic[1+is] <- 0
gam[1+is] <- alpha
pownrc[1+is] <- 0
powc[1+is] <- alpha
probs <- exp(m*q1l + n*q2l)
pownr <- 0
pow <- powc[1+0]*probs

for (is in 1:(nn-1))                               
{  ixl <- max(0,is-n)
ixu <- min(is,m)
hyp0[2+ixu] <- 0
hypa[2+ixu] <- 0
for (j in 1:(ixu-ixl+1))
{ ix <- ixu-j
hl <- fakl[1+m]-fakl[1+ix+1]-fakl[1+m-ix-1]+fakl[1+n]-fakl[1+is-ix-1]-fakl[1+n-is+ix+1]
hyp0[2+ix] <- exp(hl + (ix+1)*rho0l - oom0l)
hypa[2+ix] <- exp(hl + (ix+1)*rhoal - oomal)
hyp0[2+ix] <- hyp0[2+ix] + hyp0[2+ix+1]
hypa[2+ix] <- hypa[2+ix] + hypa[2+ix+1]  }
oom0l <- oom0l + log(hyp0[2+ixl-1])
oomal <- oomal + log(hypa[2+ixl-1])
for (ix in ixl:(ixu-1))
{ hyp0[2+ix] <- hyp0[2+ix] / hyp0[2+ixl-1]
hypa[2+ix] <- hypa[2+ix] / hypa[2+ixl-1]  }
hyp0[2+ixl-1] <- 1
hypa[2+ixl-1] <- 1
k <- ixu+1
k <- k - 1
size <- hyp0[2+k]
while(size-alpha <= 0)
{  k <- k - 1
size <- hyp0[2+k]  }
ic[1+is] <- k+1
gam[1+is] <- (alpha-hyp0[2+ic[1+is]]) / (hyp0[2+k]-hyp0[2+ic[1+is]])
pownrc[1+is] <- hypa[2+ic[1+is]]
powc[1+is] <- pownrc[1+is]*(1-gam[1+is]) + gam[1+is]*hypa[2+k]
probs <- 0
for (ix in ixl:ixu)
{ bl <- fakl[1+n]-fakl[1+is-ix]-fakl[1+n-is+ix] + 
  (is-ix)*p2l + (n-is+ix)*q2l + fakl[1+m] - fakl[1+ix] - 
  fakl[1+m-ix] + ix*p1l + (m-ix)*q1l
probs <- probs + exp(bl)   }

pownr <- pownr + pownrc[1+is]*probs
pow <- pow + powc[1+is]*probs        }  

ic[1+nn] <- n
gam[1+nn] <- alpha
pownrc[1+nn] <- 0
powc[1+nn] <- alpha
probs <- exp(m*p1l + n*p2l)
pow <- pow + powc[1+nn]*probs 

cat(" M = ",m,"   N = ",n,"   EPS = ",eps,"\n",
    "ALPHA = ",alpha,"   P1 = ",p1,"   P2 = ",p2,"\n","\n",
    "POWNR = ",pownr,"        POW = ",pow)

erg <- c(pownr,pow)

}
