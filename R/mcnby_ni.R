mcnby_ni <- function(N,DEL0,K1,K2,K3,NSUB,SW,ALPHA,MAXH)    {

PPOST <- posterior3(N,DEL0,K1,K2,K3,NSUB,C,G)
   NPL_N0<-  NPLCRIT(N,ALPHA,PPOST) 
   RES_SIZE <- FINDSIZE2(N,ALPHA, DEL0,SW,NPL_N0) 
   SIZE_UNC <- RES_SIZE[1] 

if (SIZE_UNC <= ALPHA) {
     cat(" N =",N," DEL0 =",DEL0," ALPHA =",ALPHA," K1 =",K1," K2 =",K2," K3 =",K3," NSUB =",NSUB," SW =",SW)
     stop              }

ALPHA0 <- ALPHA  
SIZE <- SIZE_UNC  
ETA_UNC <- RES_SIZE[2] 

   while(SIZE >= ALPHA)   {
   ALPHA0 <-  ALPHA0 - .01 
   NPL_N0<-  NPLCRIT(N,ALPHA0,PPOST) 
   RES_SIZE <-  FINDSIZE2(N,ALPHA0, DEL0,SW,NPL_N0) 
   SIZE <-  RES_SIZE[1] 
                         } 
ALPHA1 <- ALPHA0  
SIZE1 <- SIZE 
ALPHA2 <- ALPHA0 + .01 
IT <- 0 
   while(IT <= MAXH)  { 
   ALPHA0<- (ALPHA1+ALPHA2)/2  
   IT<- IT+1 
   NPL_N0<- NPLCRIT(N,ALPHA0,PPOST) 
   RES_SIZE <- FINDSIZE2(N,ALPHA0, DEL0,SW,NPL_N0) 
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
cat(" N =",N," DEL0 =",DEL0," ALPHA =",ALPHA," K1 =",K1," K2 =",K2," K3 =",K3," NSUB =",NSUB," SW =",SW,
    "\n", "ALPHA0 =", ALPHA0, "SIZE0 =", SIZE0, "SIZE_UNC =", SIZE_UNC)

NPL_N0<-  NPLCRIT(N,ALPHA0,PPOST) 
ETA <-  matrix(1,7,1) 
ETA[1]<- .0002  
ETA[2]<- .002  
ETA[3]<- .02  
ETA[4]<- .2 
ETA[5]<- .3  
ETA[6]<- .5  
ETA[7]<- .8 

POW <-  POW_NULLALT(N,ETA,NPL_N0) 
cat("\n"," POW =",POW)                 
#return("") 
                              }