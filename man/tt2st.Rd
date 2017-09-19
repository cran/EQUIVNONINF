\name{tt2st}
\alias{tt2st}
\title{Critical constants and power against the null alternative of the two-sample t-test for equivalence 
with an arbitrary, maybe nonsymmetric choice of the limits of the equivalence range}
\description{The function computes the critical constants defining the uniformly most powerful 
invariant test for the problem 
\eqn{(\xi-\eta)/\sigma \le -\varepsilon_1} or \eqn{(\xi-\eta)/\sigma \ge \varepsilon_2}
versus \eqn{-\varepsilon_1 < (\xi-\eta)/\sigma < \varepsilon_2}, with \eqn{\xi} and \eqn{\eta} denoting
the expected values of two normal distributions with common variance \eqn{\sigma^2} from which independent
samples are taken.
In addition, tt2st outputs the power against the null alternative \eqn{\xi = \eta}. 
}
\usage{
tt2st(m,n,alpha,eps1,eps2,tol,itmax) 
}
\arguments{
   \item{m}{size of the sample from \eqn{{\cal N}(\xi,\sigma^2)}}
   \item{n}{size of the sample from \eqn{{\cal N}(\eta,\sigma^2)}}
   \item{alpha}{significance level}
   \item{eps1}{absolute value of the lower equivalence limit to \eqn{(\xi-\eta)/\sigma} }
   \item{eps2}{upper equivalence limit to \eqn{(\xi-\eta)/\sigma} }
   \item{tol}{tolerable deviation from \eqn{\alpha} of the
              rejection probability at either boundary of
              the hypothetical equivalence interval}
   \item{itmax}{maximum number of iteration steps}          
}
\value{
   \item{m}{size of the sample from \eqn{{\cal N}(\xi,\sigma^2)}}
   \item{n}{size of the sample from \eqn{{\cal N}(\eta,\sigma^2)}}
   \item{alpha}{significance level}
   \item{eps1}{absolute value of the lower equivalence limit to \eqn{(\xi-\eta)/\sigma} }
   \item{eps2}{upper equivalence limit to \eqn{(\xi-\eta)/\sigma} }
  \item{IT}{number of iteration steps performed until reaching 
          the stopping criterion corresponding to TOL}
\item{C1}{left-hand limit of the critical interval for the two-sample \eqn{t}-statistic}
\item{C2}{right-hand limit of the critical interval for the two-sample \eqn{t}-statistic}
\item{ERR1}{deviation of the rejection probability from \eqn{\alpha}
            under \eqn{(\xi-\eta)/\sigma= -\varepsilon_1} }
\item{ERR2}{deviation of the rejection probability from \eqn{\alpha}
            under \eqn{(\xi-\eta)/\sigma= \varepsilon_2} }
\item{POW0}{power of the UMPI test against the alternative \eqn{\xi = \eta} }
}
\references{Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 6.1.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\note{If the output value of ERR2 is NA, the deviation of the rejection probability at the right-hand
      boundary of the hypothetical equivalence interval from \eqn{\alpha} is smaller than the smallest
      real number representable in R.}
\examples{
tt2st(12,12,0.05,0.50,1.00,1e-10,50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{two-sample t-test for equivalence}
\keyword{standardized difference of means}
\keyword{homoskedasticity}
\keyword{noncentral t-distribution}
\keyword{uniformly most powerful invariant test}
\keyword{power}
