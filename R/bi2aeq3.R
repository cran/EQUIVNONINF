bi2aeq3 <-function(m,n,rho1,rho2,alpha,sw,tolrd,tol,maxh)     {

fakl <- rep(NA,2000)

fakl[1+0] <- 0
nn <- m+n

for (i in 1:(nn-1))
{ qi <- i
fakl[1+i] <- fakl[1+i-1] + log(qi)  }

alph_0 <- alpha
size <- 0
nhst <- 0

cat(" m = ",m,"   n = ",n,"   rho1 = ",rho1,"   rho2 = ",rho2,
    "\n","alpha = ",alpha,"sw = ",sw,"   tolrd = ",tolrd,"   tol = ",tol,"  maxh = ",maxh,
    file="aeq3.out",append=FALSE)

while (size <= alpha)
{ alph_0 <- alph_0 + .01
size <- rejmaxaeq(m,n,alph_0,rho1,rho2,fakl,sw,tolrd,alpha)
cat("\n","\n","alph_0 =",alph_0,"   NHST =",nhst,"   SIZE =",size,
    file="aeq3.out",append=TRUE) 
                                                                            }

nhst <- 0
alph_1 <- alph_0 - .01
alph_2 <- alph_0

repeat                                              
{ nhst <- nhst + 1
  alph_0 <- (alph_1 + alph_2) / 2
  size <- rejmaxaeq(m,n,alph_0,rho1,rho2,fakl,sw,tolrd,alpha)
  cat("\n","alph_0 =",alph_0,"   NHST =",nhst,"   SIZE =",size,
      file="aeq3.out",append=TRUE)
  if (abs(size-alpha) < tol || nhst >= maxh || (size > alpha & size < alpha+.001)) break
  if (size > alpha) alph_2 <- alph_0
  if (size < alpha) alph_1 <- alph_0    }         

alph_0 <- alph_1
size <- rejmaxaeq(m,n,alph_0,rho1,rho2,fakl,sw,tolrd,alpha)
cat("\n","alph_0 =",alph_0,"   NHST =",nhst,"   SIZE =",size,
    file="aeq3.out",append=TRUE)                              

file.show("aeq3.out")

                                      }                          
