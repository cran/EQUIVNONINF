bi2ste3 <-
function(m,n,eps,alpha,sw,tolrd,tol,maxh)   {

outste3 <- tempfile()

cat(" m = ",m,"   n = ",n,"   eps = ",eps,"   alpha = ",alpha,
    "\n","sw = ",sw,"   tolrd = ",tolrd,"   tol = ",tol,"   maxh = ",maxh,"\n",
    file=outste3,append=FALSE)
  

fakl <- rep(0,2000)
hyp0 <- matrix(rep(0,2000*2000),nrow=2000)

rho0l <- log(1-eps)
nn <- m+n
fakl[1+0] <- 0

for (i in 1:2000)
{ qi <- i
  fakl[1+i] <- fakl[1+i-1] + log(qi)  }

oom0l <- 0

for (is in 1:(nn-1))                                              
{  ixl <- max(0,is-n)
  ixu <- min(is,m)
  hyp0[2+is,2+ixu] <- 0
  
  for (j in 1:(ixu-ixl+1))                                        
     { ix <- ixu - j
       hl <- fakl[1+m] - fakl[1+ix+1] - fakl[1+m-ix-1] +
             fakl[1+n] - fakl[1+is-ix-1] - fakl[1+n-is+ix+1]
       hyp0[2+is,2+ix] <- exp(hl + (ix+1)*rho0l - oom0l)
       hyp0[2+is,2+ix] <- hyp0[2+is,2+ix] + hyp0[2+is,2+ix+1]
                                                               }  
  oom0l <- oom0l + log(hyp0[2+is,2+ixl-1])
  
  for (ix in ixl:(ixu-1))                                         
  { hyp0[2+is,2+ix] <- hyp0[2+is,2+ix] / hyp0[2+is,2+ixl-1]
                                                               }  
  
  hyp0[2+is,2+ixl-1] <- 1                                       }                         #  Ende DO 4

alph_0 <- alpha
nhst <- 0
size <- rejmax(m,n,eps,sw,tolrd,alph_0,fakl,hyp0)
cat("\n","alph_0 =",alph_0,"   NHST =",nhst,"   SIZE =",size,file=outste3,append=TRUE)

while (size <= alpha)
     { alph_0 <- alph_0 + .01
       size <- rejmax(m,n,eps,sw,tolrd,alph_0,fakl,hyp0)
       cat("\n","alph_0 =",alph_0,"   NHST =",nhst,
           "   SIZE =",size,file=outste3,append=TRUE)  }

nhst <- 0
alph_1 <- alph_0 - .01
alph_2 <- alph_0

repeat                                               
{  nhst <- nhst + 1
   alph_0 <- (alph_1 + alph_2) / 2
   size <- rejmax(m,n,eps,sw,tolrd,alph_0,fakl,hyp0)
   cat("\n","alph_0 =",alph_0,"   NHST =",nhst,"   SIZE =",size,file=outste3,append=TRUE)
   if (abs(size-alpha) < tol || nhst >= maxh) break
   if (abs(size-alpha) >= tol && nhst < maxh)
     { if (size > alpha && size < (alpha + 0.001) ) break  }
   if (size > alpha) alph_2 <- alph_0
   if (size < alpha) alph_1 <- alph_0    }         

size <- rejmax(m,n,eps,sw,tolrd,alph_0,fakl,hyp0)
cat("\n","alph_0 =",alph_0,"   NHST =",nhst,"   SIZE =",size,file=outste3,append=TRUE)

file.show(outste3)

  }
