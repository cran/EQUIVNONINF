\name{tt1st}
\alias{tt1st}
\title{Critical constants and power against the null alternative of the one-sample t-test for equivalence 
with an arbitrary, maybe nonsymmetric choice of the limits of the equivalence range}
\description{The function computes the critical constants defining the uniformly most powerful 
invariant test for the problem 
\eqn{\delta/\sigma_D \le \theta_1} or \eqn{\delta/\sigma_D \ge \theta_2}
versus \eqn{\theta_1 < \delta/\sigma_D < \theta_2}, with \eqn{(\theta_1,\theta_2)} as a 
fixed nondegenerate interval on the real line. 
In addition, tt1st outputs the power against the null alternative \eqn{\delta = 0}. 
}
\usage{
tt1st(n,alpha,theta1,theta2,tol,itmax) 
}
\arguments{
   \item{n}{sample size}
   \item{alpha}{significance level}
   \item{theta1}{lower equivalence limit to \eqn{\delta/\sigma_D} }
   \item{theta2}{upper equivalence limit to \eqn{\delta/\sigma_D} }
   \item{tol}{tolerable deviation from \eqn{\alpha} of the
              rejection probability at either boundary of
              the hypothetical equivalence interval}
   \item{itmax}{maximum number of iteration steps}          
}
\value{
   \item{n}{sample size}
   \item{alpha}{significance level}
   \item{theta1}{lower equivalence limit to \eqn{\delta/\sigma_D} }
   \item{theta2}{upper equivalence limit to \eqn{\delta/\sigma_D} }
  \item{IT}{number of iteration steps performed until reaching 
          the stopping criterion corresponding to TOL}
\item{C1}{left-hand limit of the critical interval for the one-sample \eqn{t}-statistic}
\item{C2}{right-hand limit of the critical interval for the one-sample \eqn{t}-statistic}
\item{ERR1}{deviation of the rejection probability from \eqn{\alpha}
            under \eqn{\delta/\sigma_D = \theta_1}}
\item{ERR2}{deviation of the rejection probability from \eqn{\alpha}
            under \eqn{\delta/\sigma_D = \theta_2}}
\item{POW0}{power of the UMPI test against the alternative \eqn{\delta = 0}}
}
\references{Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 5.3.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\note{If the output value of ERR2 is NA, the deviation of the rejection probability at the right-hand
      boundary of the hypothetical equivalence interval from \eqn{\alpha} is smaller than the smallest
      real number representable in R.}
\examples{
tt1st(36,0.05, -0.4716,0.3853,1e-10,50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{paired t-test for equivalence}
\keyword{standardized difference of means}
\keyword{noncentral t-distribution}
\keyword{uniformly most powerful invariant test}
\keyword{power}
