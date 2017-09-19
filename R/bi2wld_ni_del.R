bi2wld_ni_del <- function(N1,N2,EPS,SW,ALPHA,MAXH)    {

TEPS <- teststat(N1,N2,EPS)
   KX_Y<-  XCRIT(N1,N2,ALPHA,TEPS) 
   RES_SIZE <- FINDSIZE(N1,N2,ALPHA, EPS,SW,KX_Y) 
   SIZE_UNC <- RES_SIZE[1] 

if (SIZE_UNC <= ALPHA) {
     cat(" N1 =",N1," N2 =",N2," EPS =",EPS," ALPHA =",ALPHA," SW =",SW," SIZE_UNC =",SIZE_UNC)
     stop              }

ALPHA0 <- ALPHA  
SIZE <- SIZE_UNC  
P2_UNC <- RES_SIZE[2] 

   while(SIZE >= ALPHA)   {
   ALPHA0 <-  ALPHA0 - .01 
   KX_Y<- XCRIT(N1,N2,ALPHA0,TEPS) 
   RES_SIZE <-  FINDSIZE(N1,N2,ALPHA0, EPS,SW,KX_Y) 
   SIZE <-  RES_SIZE[1] 
                         } 
ALPHA1 <- ALPHA0  
SIZE1 <- SIZE 
ALPHA2 <- ALPHA0 + .01 
IT <- 0 
   while(IT <= MAXH)  { 
   ALPHA0<- (ALPHA1+ALPHA2)/2  
   IT<- IT+1 
   KX_Y<- XCRIT(N1,N2,ALPHA0,TEPS) 
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
ALPHA0<- ALPHA1  
SIZE0<- SIZE1 

KX_Y_final<- XCRIT(N1,N2,ALPHA0,TEPS)
U_ALC <- qnorm(1-ALPHA0)
cat("   KX_Y_final",KX_Y_final)
  ERR_IND <- 0
     for(Y_ in 1:(N2+1))  {
          if (KX_Y_final[Y_] <= N1) {
             for(X_ in (KX_Y_final[Y_]+1):(N1+1)) { 
             if (TEPS[X_,Y_]  <= U_ALC) {ERR_IND <- 1 }
                                                  } 
                                    }
                           }
cat(" N1 =",N1," N2 =",N2," EPS =",EPS," ALPHA =",ALPHA," SW =",SW,
 "\n", "ALPHA0 =", ALPHA0, "SIZE0 =", SIZE0, "SIZE_UNC =", SIZE_UNC, "ERR_IND =", ERR_IND)

#return("") 
                              }


