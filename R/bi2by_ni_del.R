bi2by_ni_del <- function(N1,N2,EPS,SW,NSUB,ALPHA,MAXH)    {

PPOST <- posterior(N1,N2,EPS,NSUB,C,G)
  KX_Y <-  XCRIT2(N1,N2,ALPHA,PPOST) 
   RES_SIZE <- FINDSIZE(N1,N2,ALPHA, EPS,SW,KX_Y) 
   SIZE_UNC <- RES_SIZE[1] 

if (SIZE_UNC <= ALPHA) {
     cat(" N1 =",N1," N2 =",N2," EPS =",EPS," ALPHA =",ALPHA," NSUB =",NSUB," SW =",SW," SIZE_UNC =",SIZE_UNC)
     stop              }

ALPHA0 <- ALPHA  
SIZE <- SIZE_UNC  
P2_UNC <- RES_SIZE[2] 

   while(SIZE >= ALPHA)   {
   ALPHA0 <-  ALPHA0 - .01 
   KX_Y <- XCRIT2(N1,N2,ALPHA0,PPOST) 
   RES_SIZE <-  FINDSIZE(N1,N2,ALPHA0, EPS,SW,KX_Y) 
   SIZE <-  RES_SIZE[1] 
                         } 
ALPHA1 <- ALPHA0  
SIZE1 <- SIZE 
ALPHA2 <- ALPHA0 + .01 
IT <- 0 
   while(IT <= MAXH)  { 
   ALPHA0 <- (ALPHA1+ALPHA2)/2  
   IT <- IT+1 
   KX_Y <- XCRIT2(N1,N2,ALPHA0,PPOST) 
   RES_SIZE <- FINDSIZE(N1,N2,ALPHA0, EPS,SW,KX_Y) 
   SIZE <- RES_SIZE[1] 
        if(SIZE < ALPHA)  { 
                ALPHA1 <- ALPHA0  
                SIZE1 <- SIZE 
                          }
        else  {
        ALPHA2 <-  ALPHA0 
                           }
                      } 

ALPHA0 <- ALPHA1  
SIZE0 <- SIZE1 
cat(" N1 =",N1," N2 =",N2," EPS =",EPS," ALPHA =",ALPHA," NSUB =",NSUB," SW =",SW,
    "\n", "ALPHA0 =", ALPHA0, "SIZE0 =", SIZE0, "SIZE_UNC =", SIZE_UNC)

#return("") 
                              }
 


