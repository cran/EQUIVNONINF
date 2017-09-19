mcnby_ni_pp <- function(N,DEL0,N10,N01) {

ppost<- ppost_singleobs(N,DEL0,N10,N01) 

cat(" n =",N,"  del0 =",DEL0,"  n10 =",N10,"  n01 =",N01,
   "\n", "PPOST =", ppost)
#return("")
                                           }