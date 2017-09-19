cf_reh_midp <- function(X1,X2,X3,alpha,SW,TOL,ITMAX)        {

cat("\n"," X1=",X1,"   X2=",X2,"   X3=", X3,"   alpha=",alpha,
"\n"," SW=",SW, " TOL=",TOL," ITMAX=",ITMAX)
ConfAsy <- asympt_confbs(X1,X2,X3,alpha)                   

thet0 <-  4*ConfAsy[1,1]^2
exactres_l <-  exact_confb_l(X1,X2,X3,alpha, SW,TOL,ITMAX, thet0)
C_l_midp <- exactres_l[1,1]
prob_ghe <- exactres_l[1,2]
IT <- exactres_l[1,3]
cat(" C_l_midp", C_l_midp)

thet0 <- 4*ConfAsy[1,3]^2
exactres_r <- exact_confb_r(X1,X2,X3,alpha,SW,TOL,ITMAX, thet0)
C_r_midp <- exactres_r[1,1]
prob_lhe <- exactres_r[1,2]
IT <- exactres_r[1,3]
cat("\n","\n"," C_l_midp", C_l_midp,"  C_r_midp", C_r_midp)  }