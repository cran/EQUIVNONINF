\name{exp1st}
\alias{exp1st}
\title{Critical constants and power against the null alternative of the    
       UMP test for equivalence of the hazard rate of a single exponential  
       distribution to some given reference value}
\description{The function computes the critical constants defining the uniformly most powerful test for the problem 
\eqn{\sigma \le 1/(1 + \varepsilon)} or \eqn{\sigma\ge (1 + \varepsilon)}
versus \eqn{1/(1 + \varepsilon) < \sigma < (1 + \varepsilon)}, 
with \eqn{\sigma} denoting the scale parameter [\eqn{\equiv} reciprocal hazard rate] of an exponential distribution. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
exp1st(alpha,tol,itmax,n,eps) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{alpha}{significance level}
   \item{tol}{tolerable deviation from \eqn{\alpha} of the
              rejection probability at either boundary of
              the hypothetical equivalence interval}
   \item{itmax}{maximum number of iteration steps}           
   \item{n}{sample size}
   \item{eps}{margin determining the hypothetical equivalence range
              symmetrically on the log-scale}
}
\value{
   \item{alpha}{significance level}
   \item{tol}{tolerable deviation from \eqn{\alpha} of the
              rejection probability at either boundary of
              the hypothetical equivalence interval}
   \item{itmax}{maximum number of iteration steps}           
   \item{n}{sample size}
   \item{eps}{margin determining the hypothetical equivalence range
              symmetrically on the log-scale}
\item{IT}{number of iteration steps performed until reaching 
          the stopping criterion corresponding to TOL}
\item{C1}{left-hand limit of the critical interval for
          \eqn{T =\sum_{i=1}^n X_i}}
\item{C2}{right-hand limit of the critical interval for
          \eqn{T =\sum_{i=1}^n X_i}}
\item{ERR1}{deviation of the rejection probability from \eqn{\alpha}
            under \eqn{\sigma = 1/(1 + \varepsilon)} }
\item{POW0}{power of the randomized UMP test against the
            alternative \eqn{\sigma = 1}}
}
\references{Wellek S: Testing statistical hypotheses of equivalence and noninferiority. 
Second edition. Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 4.2.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
exp1st(0.05,1.0e-10,100,80,0.3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{exponential distribution}
\keyword{one-sample problem}
\keyword{UMP test for equivalence}
\keyword{power}
