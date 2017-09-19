pwexact <-
  function(m,n,eps,alpha,p1,p2)  {
    
    rho0l <- log(1-eps)
    p1l <- log(p1)
    q1l <- log(1-p1)
    p2l <- log(p2)
    q2l <- log(1-p2)
    rhoal <- p1l + q2l - q1l - p2l  
    
    ic <- rep(NA,2000)
    gam <- rep(NA,2000)
    powc <- rep(NA,2000)
    hyp0 <- rep(NA,2000)
    hypa <- rep(NA,2000)
    fakl <- rep(NA,2000)
    
    fakl[1+0] <- 0
    for (i in 1:2000)
    { qi <- i
    fakl[1+i] <- fakl[1+i-1] + log(qi) }
    
    nn <- m + n
    i <- 0
    oom0l <- 0
    oomal <- 0
    ic[1+i] <- 0
    gam[1+i] <- alpha
    powc[1+i] <- alpha
    probs <- exp(m*q1l + n*q2l)
    pow <- powc[1+0]*probs
    
    for (i in 1:(nn-1))                               
    {  ixl <- max(0,i-n)
    ixu <- min(i,m)
    hyp0[2+ixu] <- 0
    hypa[2+ixu] <- 0
    for (j in 1:(ixu-ixl+1))
    { ix <- ixu-j
    hl <- fakl[1+m]-fakl[1+ix+1]-fakl[1+m-ix-1]+fakl[1+n]-fakl[1+i-ix-1]-fakl[1+n-i+ix+1]
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
    ic[1+i] <- k+1
    gam[1+i] <- (alpha-hyp0[2+ic[1+i]]) / (hyp0[2+k]-hyp0[2+ic[1+i]])
    powc[1+i] <- hypa[2+ic[1+i]]*(1-gam[1+i]) + gam[1+i]*hypa[2+k]
    probs <- 0
    for (ix in ixl:ixu)
    { bl <- fakl[1+n]-fakl[1+i-ix]-fakl[1+n-i+ix] + 
      (i-ix)*p2l + (n-i+ix)*q2l + fakl[1+m] - fakl[1+ix] - 
      fakl[1+m-ix] + ix*p1l + (m-ix)*q1l
    probs <- probs + exp(bl)   }
    
    pow <- pow + powc[1+i]*probs        }   
    
    ic[1+nn] <- n
    gam[1+nn] <- alpha
    powc[1+nn] <- alpha
    probs <- exp(m*p1l + n*p2l)
    pwexact <- pow + powc[1+nn]*probs      }

rejmax <-
function(m,n,eps,sw,tolrd,alph_0,fakl,hyp0)   {
  
  ic <- rep(NA,2000) 
  rejpbc <- rep(NA,2000)
  
  nn <- m+n
  
  for (is in 1:(nn-1))                        
     {  ixl <- max(0,is-n)
        ixu <- min(is,m)
        k <- ixu + 1
        k <- k - 1
        size <- hyp0[2+is,2+k]
        while(size-alph_0 <= 0)
           {  k <- k - 1
              size <- hyp0[2+is,2+k]  }
        
        ic[1+is] <- k + 1                  }    
  
  itmax <- 1/sw - 1
  rejmax <- 0
  p2 <- tolrd
  p1 <- (1-eps)*p2/(1-p2)
  p1 <- p1/(1+p1)
  p2_0 <- p2
  p1_0 <- p1
  pbrj <- rjpbnc7(m,n,p2,p1,ic,fakl,hyp0)
  rejmax <- max(rejmax,pbrj)
  
  for (it in 1:itmax)                            
  { p2 <- it*sw
    p1rd <- (1-eps)*p2 / (1-p2)
    p1rd <- p1rd / (1+p1rd)
    p1 <- p1rd
    pbrj <- rjpbnc7(m,n,p2,p1,ic,fakl,hyp0)
    rejmax <- max(rejmax,pbrj)
                                                  }    
  
  p2 <- 1 - tolrd
  p1 <- (1-eps) * p2 / (1-p2)
  p1 <- p1 / (1+p1)
  pbrj <- rjpbnc7(m,n,p2,p1,ic,fakl,hyp0)
  rejmax <- max(rejmax,pbrj)     
  
  
                                                }

rjpbnc7 <-
function(m,n,p2,p1,ic,fakl,hyp0)    {
  
  rejpb_c7 <- rep(NA,2000)
  
  nn <- m + n
  p1l <- log(p1)
  q1l <- log(1-p1)
  p2l <- log(p2)
  q2l <- log(1-p2)
  
  rjpbnc7 <- 0
  
  for (is in 1:(nn-1))                            
     {  ixl <- max(0,is-n)
        ixu <- min(is,m)
        rejpb_c7[1+is] <- hyp0[2+is,2+ic[1+is]]
        probs <- 0
        for (ix in ixl:ixu)
           { bl <- fakl[1+n]-fakl[1+is-ix]-fakl[1+n-is+ix] + 
                   (is-ix)*p2l + (n-is+ix)*q2l + fakl[1+m] - fakl[1+ix] - 
                   fakl[1+m-ix] + ix*p1l + (m-ix)*q1l
             probs <- probs + exp(bl)   }
  
  rjpbnc7 <- rjpbnc7 + rejpb_c7[1+is]*probs  
                                                 }    
  return(rjpbnc7)
}

teststat <- function(N1,N2,EPS)  {     
testeps <- matrix(1,N1+1,N2+1)
testeps[1,1] <- 0
testeps[N1+1,N2+1] <- 10

Y<-0
    for(X in 1:(N1-1))   {
    T<-((X/N1-Y/N2)+EPS)/sqrt((1/N1)*(X/N1)*(1-X/N1)+
                             (1/N2)*(Y/N2)*(1-Y/N2))
    testeps[X+1,Y+1] <- T }
    
    testeps[N1+1,Y+1] <- 10
for(Y in 1:(N2-1))   {
    for(X in 0:N1)   {
    T<-((X/N1-Y/N2)+EPS)/sqrt((1/N1)*(X/N1)*(1-X/N1)+
                             (1/N2)*(Y/N2)*(1-Y/N2))
    testeps[X+1,Y+1] <- T    }
                    }
Y<-N2
    testeps[1,Y+1] <- -10
    for(X in 1:(N1-1))   {
    T<-((X/N1-Y/N2)+EPS)/sqrt((1/N1)*(X/N1)*(1-X/N1)+
                             (1/N2)*(Y/N2)*(1-Y/N2))
    testeps[X+1,Y+1] <- T
    }
return(testeps)
                                 }
								 
XCRIT <- function(N1,N2,ALPHA,TEPS)    {
U_ALC <- qnorm(1-ALPHA)
KX_Y <- matrix(1,N2+1,1)
for(Y_ in 1:(N2+1))    {
   X_ <- 1
   while(TEPS[X_,Y_] <  U_ALC & X_ <= N1) {
   X_ <- X_+1                     }
KX_Y[Y_] <- X_-1  
if (KX_Y[Y_] == N1 & TEPS[N1+1,Y_] > U_ALC ){
KX_Y[Y_] <- N1+1}
}
return(KX_Y)                        }

FINDSIZE <- function(N1,N2,ALPHA,EPS,SW,KX_Y)      {

PBY <- matrix(1,N2+1,1)

SIZE<-0  
SIZE_ <- 0  
PI2_ <- 0

pi2<- EPS

while (pi2 <= 1 -SW) { 
pi2<- pi2+SW
   PBY[1]=pbinom(0,N2,pi2)
   for(Y_ in 2:(N2+1)){ 
   Y<- Y_-1    
   PBY[Y_]<- pbinom(Y,N2,pi2) - pbinom(Y-1,N2,pi2)  
                      }
   pi1<- pi2-EPS
   PROBRJ<-0 
   for(Y_ in 1:(N2+1))    { 
     Y<- Y_-1
        if (KX_Y[Y_]==0)   {
                PBXGEK_Y<-1   } 

        if (KX_Y[Y_]== N1+1)   {
                PBXGEK_Y<-0   }
        if (KX_Y[Y_] >= 1 & KX_Y[Y_] <= N1) {
                PBXGEK_Y<- 1-pbinom(KX_Y[Y_]-1,N1,pi1)} 
                
           PBYEQY<- PBY[Y_] 
             if (min(PBXGEK_Y,PBYEQY) <= 0) { 
                        INCR = 0  }
             if (min(PBXGEK_Y,PBYEQY) > 0) {
                        LPBX<- log(PBXGEK_Y)  
                        LPBY<- log(PBYEQY) 
                        INCR<- exp(LPBX+LPBY) }                                     
     PROBRJ<- PROBRJ+INCR }
     SIZE<-max(SIZE_,PROBRJ) 
        if (SIZE > SIZE_) {  
            SIZE_ <- SIZE  
            PI2_ <- pi2 }
                                 }
                                     
RESULTS_SIZE <- matrix(1,2,1)  
RESULTS_SIZE[1] <- SIZE   
RESULTS_SIZE[2] <- PI2_ 
return(RESULTS_SIZE)                            
                        }

asympt_confbs<- function(X1,X2,X3,alpha)   {
n<- X1+X2+X3
ConfAsy <-  matrix(1,1,3)
pi1h <-  X1/n 
pi2h <-  X2/n 
pi3h <-  X3/n
thet_h <-  pi2h**2/(pi1h*pi3h)
stderr<-  sqrt((1/n)*((1-pi2h)/(pi1h*pi3h) + 4/pi2h))
om_h <-  sqrt(thet_h)/2
C_l <-  log(thet_h) - qnorm(1-alpha)*stderr
C_r <-  log(thet_h) + qnorm(1-alpha)*stderr
C_l_scaled <-  exp(.5*(C_l - log(4)))
C_r_scaled <-  exp(.5*(C_r - log(4)))

ConfAsy[1,1] <-  C_l_scaled 
ConfAsy[1,2] <-  om_h
ConfAsy[1,3] <-  C_r_scaled
return(ConfAsy)
                                          }

ex_prob_ge<- function(X1,X2,X3,thet)    {
S <-  2*X1 + X2 
N<- X1+X2+X3
NB<- floor(S/2)-max(0,S-N)+1

B <- matrix(1,1,NB)  
PRB <-  matrix(1,1,NB)

for(K in 1:NB) { B[1,K] <- S-2*floor(S/2)+2*(K-1) }
K<- 1 
while(B[1,K] <= X2) { K<- K+1 } 
KX2<-  K-1
        K <-  1
        CL<- lgamma(N+1)-lgamma((S-B[1,K])/2+1) -lgamma(B[1,K]+1) -
           lgamma(N-B[1,K]/2- S/2+1)
        ARGEXP_U <-  CL+(B[1,K]/2)*log(thet)

        for(K in 2:NB)   {
        CL<- lgamma(N+1)-lgamma((S-B[1,K])/2+1) -lgamma(B[1,K]+1) -
           lgamma(N-B[1,K]/2- S/2+1)
        ARGEXP_U <-  max(ARGEXP_U,  CL+(B[1,K]/2)*log(thet))
                          }

        SHIFTL <-  min(0,700-ARGEXP_U)

        for(K in 1:NB)   {
        CL<- lgamma(N+1)-lgamma((S-B[1,K])/2+1) -lgamma(B[1,K]+1) -
           lgamma(N-B[1,K]/2- S/2+1)
        PRB[1,K]<-  exp(CL+(B[1,K]/2)*log(thet) + SHIFTL)
                         }

      for(K in 2:NB)     {
      PRB[1,K]<-  PRB[1,K]+PRB[1,K-1]
                         }

      for(K in 1:NB)   {
      PRB[1,K]<- PRB[1,K]/PRB[1,NB]
                       }

prob_ge <-  1- PRB[1,KX2-1]
return(prob_ge)                         } 

ex_prob_le<- function(X1,X2,X3,thet)    {
S <-  2*X1 + X2 
N<- X1+X2+X3
NB<- floor(S/2)-max(0,S-N)+1

B <-  matrix(1,1,NB)  
PRB <-  matrix(1,1,NB)

for(K in 1:NB) { B[1,K] <- S-2*floor(S/2)+2*(K-1) }

K<- 1 
while(B[1,K] <= X2) { K<- K+1 } 
KX2<-  K-1
        K <-  1
        CL<- lgamma(N+1)-lgamma((S-B[1,K])/2+1) -lgamma(B[1,K]+1) -
           lgamma(N-B[1,K]/2- S/2+1)
        ARGEXP_U <-  CL+(B[1,K]/2)*log(thet)

        for(K in 2:NB)   {
        CL<- lgamma(N+1)-lgamma((S-B[1,K])/2+1) -lgamma(B[1,K]+1) -
           lgamma(N-B[1,K]/2- S/2+1)
        ARGEXP_U <-  max(ARGEXP_U,  CL+(B[1,K]/2)*log(thet))
                          }
        SHIFTL <-  min(0,700-ARGEXP_U)

               for(K in 1:NB)   {
        CL<- lgamma(N+1)-lgamma((S-B[1,K])/2+1) -lgamma(B[1,K]+1) -
           lgamma(N-B[1,K]/2- S/2+1)
        PRB[1,K]<-  exp(CL+(B[1,K]/2)*log(thet) + SHIFTL)
                                }

      for(K in 2:NB)   {
      PRB[1,K]<-  PRB[1,K]+PRB[1,K-1]
                       }

      for(K in 1:NB)   {
      PRB[1,K]<- PRB[1,K]/PRB[1,NB]
                       }

prob_le <-  PRB[1,KX2]
return(prob_le)
                                       }
									   
exact_confb_l <- function(X1,X2,X3,alpha, SW, TOL, ITMAX, thet0)    {
thet<- thet0  
IT <- 0
prob_ge <- ex_prob_ge(X1,X2,X3, thet)

if (prob_ge < alpha)           {
    while (prob_ge < alpha)                { 
    thet <-  thet + SW
    prob_ge <-  ex_prob_ge(X1,X2,X3, thet) }
    thet1 <-  thet-SW 
    thet2 <-  thet  
    while (IT < ITMAX)     { 
    	thet <-  (thet1+thet2)/2 
   	prob_ge <-  ex_prob_ge(X1,X2,X3, thet)
       		if (prob_ge < alpha -TOL) { 
           	thet1 <-  thet
                IT <-  IT+1 }
       		if (prob_ge > alpha +TOL) { 
           	thet2 <-  thet
                IT <-  IT+1 } 
                else  {IT <-  IT+1}                              
                                            } 
                                }

if (prob_ge > alpha)           {
    while (prob_ge > alpha)                { 
    thet <-  thet - SW
    prob_ge <-  ex_prob_ge(X1,X2,X3, thet) }
    thet1 <-  thet 
    thet2 <-  thet+SW  
    while (IT < ITMAX)     { 
    	thet <-  (thet1+thet2)/2 
   	prob_ge <-  ex_prob_ge(X1,X2,X3, thet)
       		if (prob_ge < alpha -TOL) { 
           	thet1 <-  thet
                IT <-  IT+1 }
       		if (prob_ge > alpha +TOL) { 
           	thet2 <-  thet
                IT <-  IT+1 } 
                else  {IT <-  IT+1}                              
                                            } 
                                }

if (prob_ge == alpha)           {thet<- thet0}

C_l_exact <- sqrt(thet)/2   
prob_ge <- ex_prob_ge(X1,X2,X3, thet)

exactres_l <-  matrix(1,1,3)
exactres_l[1,1] <-  C_l_exact
exactres_l[1,2] <-  prob_ge
exactres_l[1,3] <-  IT

return(exactres_l)                                }

exact_confb_r <- function(X1,X2,X3,alpha, SW, TOL, ITMAX, thet0)    {
thet<- thet0 
IT <-  0
prob_le <-  ex_prob_le(X1,X2,X3,thet)

if (prob_le < alpha)           {
    while (prob_le < alpha)                { 
    thet <-  thet - SW
    prob_le <-  ex_prob_le(X1,X2,X3, thet) }
    thet1 <-  thet 
    thet2 <-  thet+SW  
    while (IT < ITMAX)     { 
    	thet <-  (thet1+thet2)/2 
   	prob_le <-  ex_prob_le(X1,X2,X3, thet)
       		if (prob_le < alpha -TOL) { 
           	thet2 <-  thet
                IT <-  IT+1 }
       		if (prob_le > alpha +TOL) { 
           	thet1 <-  thet
                IT <-  IT+1 } 
                else  {IT <-  IT+1}                              
                                            } 
                                }

if (prob_le > alpha)           {
    while (prob_le > alpha)                { 
    thet <-  thet + SW
    prob_le <-  ex_prob_le(X1,X2,X3, thet) }
    thet1 <-  thet-SW 
    thet2 <-  thet  
    while (IT < ITMAX)     { 
    	thet <-  (thet1+thet2)/2 
   	prob_le <-  ex_prob_le(X1,X2,X3, thet)
       		if (prob_le < alpha -TOL) { 
           	thet2 <-  thet
                IT <-  IT+1 }
       		if (prob_le > alpha +TOL) { 
           	thet1 <-  thet
                IT <-  IT+1 } 
                else  {IT <-  IT+1}                              
                                            } 
                                }
if (prob_le == alpha)    {thet<- thet0}
C_r_exact <-  sqrt(thet)/2   
prob_le <-  ex_prob_le(X1,X2,X3, thet)

exactres_r <-  matrix(1,1,3)
exactres_r[1,1] <-  C_r_exact
exactres_r[1,2] <-  prob_le
exactres_r[1,3] <-  IT

return(exactres_r)                                }

ex_prob_ghe<- function(X1,X2,X3,thet)       {
S <-  2*X1 + X2 
N<- X1+X2+X3
NB<- floor(S/2)-max(0,S-N)+1

B <- matrix(1,1,NB)  
PRB <-  matrix(1,1,NB)

for(K in 1:NB) { B[1,K] <- S-2*floor(S/2)+2*(K-1) }
K<- 1 
while(B[1,K] <= X2) { K<- K+1 } 
KX2<-  K-1
        K <-  1
        CL<- lgamma(N+1)-lgamma((S-B[1,K])/2+1) -lgamma(B[1,K]+1) -
           lgamma(N-B[1,K]/2- S/2+1)
        ARGEXP_U <-  CL+(B[1,K]/2)*log(thet)

        for(K in 2:NB)   {
        CL<- lgamma(N+1)-lgamma((S-B[1,K])/2+1) -lgamma(B[1,K]+1) -
           lgamma(N-B[1,K]/2- S/2+1)
        ARGEXP_U <-  max(ARGEXP_U,  CL+(B[1,K]/2)*log(thet))
                          }

        SHIFTL <-  min(0,700-ARGEXP_U)

        for(K in 1:NB)   {
        CL<- lgamma(N+1)-lgamma((S-B[1,K])/2+1) -lgamma(B[1,K]+1) -
           lgamma(N-B[1,K]/2- S/2+1)
        PRB[1,K]<-  exp(CL+(B[1,K]/2)*log(thet) + SHIFTL)
                         }

      for(K in 2:NB)     {
      PRB[1,K]<-  PRB[1,K]+PRB[1,K-1]
                         }

      for(K in 1:NB)   {
      PRB[1,K]<- PRB[1,K]/PRB[1,NB]
                       }

prob_ghe <-  1- (PRB[1,KX2-1]+PRB[1,KX2])/2
return(prob_ghe)                              } 

ex_prob_lhe<- function(X1,X2,X3,thet)        {
S <-  2*X1 + X2 
N<- X1+X2+X3
NB<- floor(S/2)-max(0,S-N)+1

B <-  matrix(1,1,NB)  
PRB <-  matrix(1,1,NB)

for(K in 1:NB) { B[1,K] <- S-2*floor(S/2)+2*(K-1) }

K<- 1 
while(B[1,K] <= X2) { K<- K+1 } 
KX2<-  K-1
        K <-  1
        CL<- lgamma(N+1)-lgamma((S-B[1,K])/2+1) -lgamma(B[1,K]+1) -
           lgamma(N-B[1,K]/2- S/2+1)
        ARGEXP_U <-  CL+(B[1,K]/2)*log(thet)

        for(K in 2:NB)   {
        CL<- lgamma(N+1)-lgamma((S-B[1,K])/2+1) -lgamma(B[1,K]+1) -
           lgamma(N-B[1,K]/2- S/2+1)
        ARGEXP_U <-  max(ARGEXP_U,  CL+(B[1,K]/2)*log(thet))
                          }
        SHIFTL <-  min(0,700-ARGEXP_U)

               for(K in 1:NB)   {
        CL<- lgamma(N+1)-lgamma((S-B[1,K])/2+1) -lgamma(B[1,K]+1) -
           lgamma(N-B[1,K]/2- S/2+1)
        PRB[1,K]<-  exp(CL+(B[1,K]/2)*log(thet) + SHIFTL)
                                }

      for(K in 2:NB)   {
      PRB[1,K]<-  PRB[1,K]+PRB[1,K-1]
                       }

      for(K in 1:NB)   {
      PRB[1,K]<- PRB[1,K]/PRB[1,NB]
                       }

prob_lhe <- (PRB[1,KX2-1]+PRB[1,KX2])/2
return(prob_lhe)                                  }

posterior <- function(N1,N2,EPS,NSUB,C,G)    {
ppost<- matrix(1,N1+1,N2+1)
A<- EPS 
B<- 1

for(X in 0:N1)                     {
   for(Y in 0:N2)              {

PROBPOST<- 0
JJ <- 1  
     while (JJ <= NSUB)  { 
	
      AA <- A+(JJ-1)*(B-A)/NSUB  
      BB <- A+JJ*(B-A)/NSUB 
      M <- (AA+BB)/2  
      L <- (BB-AA)/2 
    
       for(J in 1:48)   { 
       ITGDR<-0 
       ITGDL<-0        
        PI2_R<-  M + C[J]*L
        PP_PI2_R<- 1- pbeta(PI2_R-EPS,X+.5,N1-X+.5)
      	dens <- dbeta(PI2_R,Y+.5,N2-Y+.5) 
  	     if(min(PP_PI2_R,dens) > 0)  {
         ITGDRLOG <- log(dens)+log(PP_PI2_R) 
         ITGDR<- exp(ITGDRLOG)      }
       
         PI2_L<-  M - C[J]*L
         PP_PI2_L<- 1- pbeta(PI2_L-EPS,X+.5,N1-X+.5)
         dens <- dbeta(PI2_L,Y+.5,N2-Y+.5) 
  	     if(min(PP_PI2_L,dens) > 0)  {
         ITGDLLOG <- log(dens)+log(PP_PI2_L) 
         ITGDL<- exp(ITGDLLOG)      }       
        PROBPOST<- PROBPOST + G[J]*(ITGDR+ITGDL)
 }  
  JJ<- JJ+1  }
PROBPOST<- PROBPOST*L+pbeta(EPS,Y+.5,N2-Y+.5) 
X_<- X+1 
Y_<- Y+1
ppost[X_,Y_]<- PROBPOST                                           }
                                           }
return(ppost)
                                                            }
															
XCRIT2 <- function(N1,N2,ALPHA,PPOST)    {
KX_Y <- matrix(1,N2+1,1)
for(Y_ in 1:(N2+1))    {
   X_ <- 1
   while(PPOST[X_,Y_] <= 1-ALPHA & X_ <= N1) {
   X_ <- X_+1                     }
KX_Y[Y_] <- X_-1  
if (KX_Y[Y_] == N1 & PPOST[N1+1,Y_] <= 1-ALPHA){
KX_Y[Y_] <- N1+1}
}
return(KX_Y)                        }	

posterior2 <- function(N1,N2,EPS,NSUB,C,G)    {
  
  RHO0 <- 1-EPS 
  ppost<- matrix(1,N1+1,N2+1)
  
  A<- 0 
  B<- 1
  
  for(X in 0:N1)                     {
    for(Y in 0:N2)              {
      
      PROBPOST<- 0
      JJ <- 1  
      while (JJ <= NSUB)  { 
        
        AA <- A+(JJ-1)*(B-A)/NSUB  
        BB <- A+JJ*(B-A)/NSUB 
        M <- (AA+BB)/2  
        L <- (BB-AA)/2 
        
        for(J in 1:48)   { 
          ITGDR<-0 
          ITGDL<-0        
          PI2_R<-  M + C[J]*L
          PI2_R_TR<- RHO0*PI2_R/((1-PI2_R)+RHO0*PI2_R)
          PP_PI2_R<- 1- pbeta(PI2_R_TR,X+.5,N1-X+.5)
          dens <- dbeta(PI2_R,Y+.5,N2-Y+.5) 
          if(min(PP_PI2_R,dens) > 0)  {
            ITGDRLOG <- log(dens)+log(PP_PI2_R) 
            ITGDR<- exp(ITGDRLOG)      }
          
          PI2_L<-  M - C[J]*L
          PI2_L_TR<- RHO0*PI2_L/((1-PI2_L)+RHO0*PI2_L);
          PP_PI2_L<- 1- pbeta(PI2_L_TR,X+.5,N1-X+.5)
          dens <- dbeta(PI2_L,Y+.5,N2-Y+.5) 
          if(min(PP_PI2_L,dens) > 0)  {
            ITGDLLOG <- log(dens)+log(PP_PI2_L) 
            ITGDL<- exp(ITGDLLOG)      }       
          PROBPOST<- PROBPOST + G[J]*(ITGDR+ITGDL)
        }  
        JJ<- JJ+1  }
      PROBPOST<- PROBPOST*L
      X_<- X+1 
      Y_<- Y+1
      ppost[X_,Y_]<- PROBPOST                                       }
  }
  return(ppost)
}

posterior3 <- function(N,DEL0,K1,K2,K3,NSUB,C,G)    {

ppost<- matrix(1,N+1,N+1)

A<- DEL0 
B<- (1+DEL0)/2

for(X1 in 0:N)                     {
   for(X2 in 0:(N-X1))              {

PROBPOST_ <- pbeta(DEL0,X2+K2,N-X2+K1+K3) 
PROBPOST<- 0
JJ <- 1  
     while (JJ <= NSUB)  { 
	
      AA <- A+(JJ-1)*(B-A)/NSUB  
      BB <- A+JJ*(B-A)/NSUB 
      M <- (AA+BB)/2  
      L <- (BB-AA)/2 
    
       for(J in 1:48)   { 
       ITGDR<-0 
       ITGDL<-0
       pmiR<-  M + C[J]*L
        Q1<- pmiR-DEL0 
      	IQ1 <- pbeta(Q1/(1-pmiR),X1+K1,N-X1-X2+K3) 
      	IQ2 <- 1 
      	IQ2_1 <- IQ2-IQ1 
      	dens <- dbeta(pmiR,X2+K2,N-X2+K1+K3) 
  	     if( min(IQ2_1,dens) > 0)  {
         ITGDRLOG <- log(dens)+log(IQ2_1) 
         ITGDR<- exp(ITGDRLOG)      }
       pmiL<-  M - C[J]*L     	
        Q1<- pmiL-DEL0 
        IQ1<-pbeta(Q1/(1-pmiL),X1+K1,N-X1-X2+K3) 
        IQ2<- 1 
        IQ2_1<- IQ2-IQ1 
        dens <- dbeta(pmiL,X2+K2,N-X2+K1+K3) 
 	      if( min(IQ2_1,dens) > 0)  {
         ITGDLLOG <- log(dens)+log(IQ2_1) 
         ITGDL<- exp(ITGDLLOG)      }        
        PROBPOST<- PROBPOST + G[J]*(ITGDR+ITGDL)
 }  
  JJ<- JJ+1  }
PROBPOST<- PROBPOST*L+PROBPOST_ 
Npl_ <- X1+1
N0_<- N-X1-X2+1
ppost[N0_,Npl_] <- PROBPOST
                                           }
                                           }
return(ppost)
                                                            }
															
NPLCRIT <- function(N,ALPHA,PPOST)    {
NPL_N0 <- matrix(1,N+1,1)
for(N0_ in 1:(N+1))    {
   NPL_ <- 1
   while(PPOST[N0_,NPL_] <= 1-ALPHA) {
   NPL_ <- NPL_+1                     }
NPL_N0[N0_] <- NPL_-1  }
return(NPL_N0)                        }

FINDSIZE2 <- function(N,ALPHA, DEL0,SW,NPL_N0)      {
SIZE<-0  
SIZE_ <- 0  
ETA_ <- DEL0 

PBN0 <- matrix(1,N+1,1) 

ETA<- DEL0 + SW 
while (ETA <= 1 -SW) { 
p0 <- 1-ETA 
pi_0 <- 1/2-DEL0/(2*ETA) 
   PBN0[1]<-pbinom(0,N,p0) 
   for( N0_ in 2:(N+1)){ 
   PBN0[N0_] <-pbinom(N0_-1,N,p0) - pbinom(N0_-2,N,p0)  
                         }
   PROBRJ<-0 
   for(N0_ in 1:N)    { 
        if (NPL_N0[N0_] <= N+1-N0_)   {
                if (NPL_N0[N0_] <=0){
                  PBNplGEK_N0<-1   } 
                else {
                  PBNplGEK_N0<- 1-pbinom(NPL_N0[N0_]-1,N+1-N0_,pi_0)} 
                PBN0EQN0<-PBN0[N0_] 
                     if (min(PBNplGEK_N0,PBN0EQN0) > 0) { 
                        LPBNpl<-log(PBNplGEK_N0)  
                        LPBN0<-log(PBN0EQN0) 
                        PROBRJ<-PROBRJ+exp(LPBNpl+LPBN0) 
                                                         }                                             
                                        }
                      }
   SIZE<-max(SIZE_,PROBRJ) 
   if (SIZE > SIZE_) {  
     SIZE_ <- SIZE  
     ETA_ <- ETA     } 
   ETA <- ETA + SW       }

RESULTS_SIZE <- matrix(1,2,1)  
RESULTS_SIZE[1] <- SIZE   
RESULTS_SIZE[2] <- ETA_ 
return(RESULTS_SIZE)                            
                        }
						
POW_NULLALT<- function(N,ETA,NPL_N0)            {
PBN0  <-  matrix(1,N+1,1)
k  <- nrow(ETA)    
POW  <- matrix(1,1,k) 

for(j in 1:k)    { 
   ETA_ <- ETA[j]   
   PI <- 1/2  
   p0 <-  1-ETA_ 
   PBN0[1] <- pbinom(0,N,p0) 
       for(N0_ in 2:(N+1))  { 
       PBN0[N0_]  <- pbinom(N0_-1,N,p0) - pbinom(N0_-2,N,p0)  
                            } 
   PROBRJ <- 0 
   for(N0_ in 1:N)  { 
        if (NPL_N0[N0_] <= N+1-N0_) {
             if (NPL_N0[N0_] == 0) {
                PBNplGEK_N0 <- 1  }
             else {               
                PBNplGEK_N0 <- 1-pbinom(NPL_N0[N0_]-1,N+1-N0_,PI) 
                PBN0EQN0 <- PBN0[N0_] 
                     if (min(PBNplGEK_N0,PBN0EQN0) > 0) { 
                        LPBNpl <- log(PBNplGEK_N0)  
                        LPBN0 <- log(PBN0EQN0) 
                        PROBRJ <- PROBRJ+exp(LPBNpl+LPBN0) 
                                                         }    
                                                   } 
                                                 }     

 POW[j] <- PROBRJ 
                                            }
                                              } 
return(POW)                                       }
ppost_singleobs <- function(N,DEL0,N10,N01)    {
 
A<- DEL0  
B<- (1+DEL0)/2 

X1<- N10  
X2<- N01  
K1<- .5  
K2 <- .5  
K3 <- .5  
NSUB <- 10 

PROBPOST_ <- pbeta(DEL0,X2+K2,N-X2+K1+K3) 

PROBPOST<- 0

JJ <- 1  
     while (JJ <= NSUB)  { 
	AA <- A+(JJ-1)*(B-A)/NSUB  
	BB <- A+JJ*(B-A)/NSUB 
	M <- (AA+BB)/2  
	L <- (BB-AA)/2 
    
       for(J in 1:48)   { 
       ITGDR<-0 
       ITGDL<-0
       pmiR<-  M + C[J]*L
        Q1<- pmiR-DEL0 
      	IQ1 <- pbeta(Q1/(1-pmiR),X1+K1,N-X1-X2+K3) 
      	IQ2 <- 1 
      	IQ2_1 <- IQ2-IQ1 
      	dens <- dbeta(pmiR,X2+K2,N-X2+K1+K3) 
  	     if( min(IQ2_1,dens) > 0)  {
         ITGDRLOG <- log(dens)+log(IQ2_1) 
         ITGDR<- exp(ITGDRLOG)      }
       pmiL<-  M - C[J]*L     	
        Q1<- pmiL-DEL0 
        IQ1<-pbeta(Q1/(1-pmiL),X1+K1,N-X1-X2+K3) 
        IQ2<- 1 
        IQ2_1<- IQ2-IQ1 
        dens <- dbeta(pmiL,X2+K2,N-X2+K1+K3) 
 	      if( min(IQ2_1,dens) > 0)  {
         ITGDLLOG <- log(dens)+log(IQ2_1) 
         ITGDL<- exp(ITGDLLOG)      }        
        PROBPOST<- PROBPOST + G[J]*(ITGDR+ITGDL)
 }  
  JJ<-JJ+1  }
PROBPOST<- PROBPOST*L+PROBPOST_ 
return(PROBPOST)                }
pwexacta <- function(m,n,alpha,q1l,q2l,p1l,p2l,rho1,rho2,rhoal,fakl)  {

rho1l <- log(rho1)
rho2l <- log(rho2)
  
hypa <- rep(0,2000)
hyp01 <- rep(0,2000)
hyp02 <- rep(0,2000)

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
      
      #  7
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
}    

powc <- alpha
probs <- exp(m*p1l + n*p2l)
pwexact <- pow + powc*probs

return(pwexact)
}
rejmaxaeq <- function(m,n,alpha0,rho1,rho2,fakl,sw,tolrd,alpha)   {
  
  ic1 <- rep(0,2000)
  ic2 <- rep(0,2000)
  hyp01 <- rep(0,2000)
  hyp02 <- rep(0,2000)
  rejpbc01 <- rep(0,2000)
  rejpbc02 <- rep(0,2000)
  
nn <- m + n
rho1l <- log(rho1)
rho2l <- log(rho2)
oom01l <- 0
oom02l <- 0

for (is in 1:(nn-1))                               
   { GT300 <- 0
     GT36  <- 0
     GT39  <- 0
     
     ixl <- max(0,is-n)
     ixu <- min(is,m)
     ixeins <- is * m/nn
     hyp01[2+ixu] <- 0
     hyp02[2+ixu] <- 0
     
     for (j in 1:(ixu-ixl+1))                     
     { ix <- ixu-j
       hl <- fakl[1+m]-fakl[1+ix+1]-fakl[1+m-ix-1]+fakl[1+n]-fakl[1+is-ix-1]-fakl[1+n-is+ix+1]
       hyp01[2+ix] <- exp(hl + (ix+1)*rho1l - oom01l)
       hyp02[2+ix] <- exp(hl + (ix+1)*rho2l - oom02l)
       hyp01[2+ix] <- hyp01[2+ix] + hyp01[2+ix+1]
       hyp02[2+ix] <- hyp02[2+ix] + hyp02[2+ix+1]
                                                      }      
     oom01l <- oom01l + log(hyp01[2+ixl-1])
     oom02l <- oom02l + log(hyp02[2+ixl-1])
     
     for (ix in ixl:(ixu-1))                              
     { hyp01[2+ix] <- hyp01[2+ix] / hyp01[2+ixl-1]
       hyp02[2+ix] <- hyp02[2+ix] / hyp02[2+ixl-1]
                                                    }     
     hyp01[2+ixl-1] <- 1
     hyp02[2+ixl-1] <- 1
     k <- ixeins
     
     if (2*k < is || rho1 != 1/rho2)
       GT300 <- 1
     else
     { hrho1k <- hyp01[2+k] - hyp01[2+k+1]
     if (hrho1k >= alpha0) GT36 <- 1       }
     
     if (GT36 == 0)
     {
       k1 <- min(k+7,is,m)                        
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
           
           if (max(alpha1,alpha2) - alpha0 > 0) break    
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
             
             delalph1 <- alpha0 - alpha1
             delalph2 <- alpha0 - alpha2
             exhyp11 <- hyp01[2+k1-2] - hyp01[2+k1-1]
             exhyp12 <- hyp01[2+k2] - hyp01[2+k2+1]
             exhyp21 <- hyp02[2+k1-2] - hyp02[2+k1-1]
             exhyp22 <- hyp02[2+k2] - hyp02[2+k2+1]
             det <- exhyp11*exhyp22 - exhyp12*exhyp21
             gamma1 <- (exhyp22*delalph1-exhyp12*delalph2)/det
             gamma2 <- (exhyp11*delalph2-exhyp21*delalph1)/det
             
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
             { GT39 <- 1
               break       }              }           
         
         if (GT39 == 1) break                     
       
                                                  }       
       
                                       }         
     
     if (GT36 == 1)                                      
     { c1 <- k
       c2 <- k
       gamma1 <- alpha0/hrho1k
       gamma2 <- gamma1         }
     
     ic1[1+is] <- k1-1                                     
     ic2[1+is] <- k2+1
     p01gtc1 <- hyp01[2+ic1[1+is]]
     p01gec2 <- hyp01[2+ic2[1+is]-1]
     p02gtc1 <- hyp02[2+ic1[1+is]]
     p02gec2 <- hyp02[2+ic2[1+is]-1]
     ic1eqc2 <- 0
     if (ic1[1+is] == ic2[1+is]) ic1eqc2 <- 1 
     
     rejpbc01[1+is] <- (p01gtc1-p01gec2)*(1-ic1eqc2)
     rejpbc02[1+is] <- (p02gtc1-p02gec2)*(1-ic1eqc2)
     
                                                       }      

itmax <- 1/sw - 1
rejmax <- 0
p2max <- 0
p2 <- tolrd

pbrj <- rjpbnc(m,n,p2,rho1,rho2,fakl,rejpbc01,rejpbc02)
rejmax <- max(rejmax,pbrj)

if(pbrj - rejmax >= 0) p2max <- p2

for(it in 1:itmax)
   { p2 <- it*sw                                            
     pbrj <- rjpbnc(m,n,p2,rho1,rho2,fakl,rejpbc01,rejpbc02)
     rejmax <- max(rejmax,pbrj)
     if(pbrj - rejmax >= 0) p2max <- p2
                                              }             
p2 <- 1-tolrd
pbrj <- rjpbnc(m,n,p2,rho1,rho2,fakl,rejpbc01,rejpbc02)
rejmax <- max(rejmax,pbrj)
if(pbrj - rejmax >= 0) p2max <- p2
return(rejmax)
                                                }
rjpbnc <- function(m,n,p2,rho1,rho2,fakl,rejpbc01,rejpbc02)    {


nn <- m + n
p1li <- rho1*p2/(1-p2)
p1li <- p1li/(1+p1li)
p1lil <- log(p1li)
q1lil <- log(1-p1li)
p1re <- rho2*p2/(1-p2)
p1re <- p1re/(1+p1re)
p1rel <- log(p1re)
q1rel <- log(1-p1re)
p2l <- log(p2)
q2l <- log(1-p2)
rjpbncli <- 0
rjpbncre <- 0

for (is in 1:(nn-1))                                       
   { probsli <- 0
     probsre <- 0
     ixl <- max(0,is-n)
     ixu <- min(is,m)
     for (ix in ixl:ixu)                                            
       { blil <- fakl[1+n]-fakl[1+is-ix]-fakl[1+n-is+ix] + 
                (is-ix)*p2l + (n-is+ix)*q2l + fakl[1+m] - fakl[1+ix] - 
                fakl[1+m-ix] + ix*p1lil + (m-ix)*q1lil
         probsli <- probsli + exp(blil)
         brel <- fakl[1+n]-fakl[1+is-ix]-fakl[1+n-is+ix] + 
                (is-ix)*p2l + (n-is+ix)*q2l + fakl[1+m] - fakl[1+ix] - 
                fakl[1+m-ix] + ix*p1rel + (m-ix)*q1rel
         probsre <- probsre + exp(brel)
                                                         }         
     rjpbncli <- rjpbncli + rejpbc01[1+is]*probsli
     rjpbncre <- rjpbncre + rejpbc02[1+is]*probsre
 
                                                      }    
rjpbnc <- max(rjpbncli,rjpbncre)
return(rjpbnc)
                                                      }