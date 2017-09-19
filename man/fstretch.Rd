\name{fstretch}
\alias{fstretch}
\title{Critical constants and power of the UMPI (uniformly most powerful invariant) test for 
dispersion equivalence of two Gaussian distributions}
\description{The function computes the critical constants defining the optimal test for the problem 
\eqn{\sigma^2/\tau^2 \le \varrho_1} or \eqn{\sigma^2/\tau^2 \ge \varrho_2}
versus \eqn{\varrho_1 < \sigma^2/\tau^2 < \varrho_2}, 
with \eqn{(\varrho_1,\varrho_2)} as a fixed nonempty interval around unity. 
}
\usage{
fstretch(alpha,tol,itmax,ny1,ny2,rho1,rho2) 
}
\arguments{
   \item{alpha}{significance level}
   \item{tol}{tolerable deviation from \eqn{\alpha} of the
              rejection probability at either boundary of
              the hypothetical equivalence interval}
   \item{itmax}{maximum number of iteration steps}           
   \item{ny1}{number of degrees of freedom of the estimator of 
              \eqn{\sigma^2}}
   \item{ny2}{number of degrees of freedom of the estimator of 
              \eqn{\tau^2}}
   \item{rho1}{lower equivalence limit to \eqn{\sigma^2/\tau^2}}
   \item{rho2}{upper equivalence limit to \eqn{\sigma^2/\tau^2}}
}
\value{
   \item{alpha}{significance level}
   \item{tol}{tolerable deviation from \eqn{\alpha} of the
              rejection probability at either boundary of
              the hypothetical equivalence interval}
   \item{itmax}{maximum number of iteration steps}           
   \item{ny1}{number of degrees of freedom of the estimator of 
              \eqn{\sigma^2}}
   \item{ny2}{number of degrees of freedom of the estimator of 
              \eqn{\tau^2}}
   \item{rho1}{lower equivalence limit to \eqn{\sigma^2/\tau^2}}
   \item{rho2}{upper equivalence limit to \eqn{\sigma^2/\tau^2}}
\item{IT}{number of iteration steps performed until reaching 
          the stopping criterion corresponding to TOL}
\item{C1}{left-hand limit of the critical interval for
        \deqn{T = \frac{n-1}{m-1} \sum_{i=1}^m (X_i-\overline{X})^2 / \sum_{j=1}^{n-1} (Y_j-\overline{Y})^2}}
\item{C2}{right-hand limit of the critical interval for
        \deqn{T = \frac{n-1}{m-1} \sum_{i=1}^m (X_i-\overline{X})^2 / \sum_{j=1}^{n-1} (Y_j-\overline{Y})^2}}
\item{ERR}{deviation of the rejection probability from \eqn{\alpha}
            under \eqn{\sigma^2/\tau^2 = \varrho_1} }
\item{POW0}{power of the UMPI test against the
            alternative \eqn{\sigma^2/\tau^2 = 1}}
}
\references{Wellek S: Testing statistical hypotheses of equivalence and noninferiority. 
            Second edition. Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 6.5.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\note{If the two independent samples under analysis are from exponential rather than Gaussian distributions, the critical constants computed by
means of fstretch with \eqn{\nu_1 = 2m}, \eqn{\nu_2 = 2n}, can be used
for testing for equivalence with respect to the ratio of hazard rates. The only difference is that the ratio of sample means rather than variances has
to be used as the test statistic then.
}
\examples{
fstretch(0.05, 1.0e-10, 50,40,45,0.5625,1.7689)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dispersion equivalence}
\keyword{two-sample problem}
\keyword{Gaussian data}
\keyword{uniformly most powerful invariant test}
\keyword{power}
