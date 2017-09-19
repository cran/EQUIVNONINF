\name{postmys}
\alias{postmys}
\title{Bayesian posterior probability of the alternative hypothesis 
in the setting of the one-sample t-test for equivalence}
\description{Evaluation of the integral appearing on the right-hand side of equation (3.6)
on p. 38 of Wellek S (2010) Testing statistical hypotheses of equivalence and 
noninferiority. Second edition}
\usage{
postmys(n,dq,sd,eps1,eps2,tol)
}
\arguments{
  \item{n}{sample size}
  \item{dq}{mean within-pair difference observed in the sample under analysis} 
  \item{sd}{square root of the sample variance of the  within-pair differences} 
  \item{eps1}{absolute value of the left-hand limit of the hypothetical equivalence range 
              for \eqn{\delta/\sigma_D} }
  \item{eps2}{right-hand limit of the hypothetical equivalence range for \eqn{\delta/\sigma_D} }
  \item{tol}{tolerance for the error induced through truncating the range of integration on the right}
}
\details{The program uses 96-point Gauss-Legendre quadrature.}
\value{
  \item{n}{sample size}
  \item{dq}{mean within-pair difference observed in the sample under analysis} 
  \item{sd}{square root of the sample variance of the  within-pair differences} 
  \item{eps1}{absolute value of the left-hand limit of the hypothetical equivalence range 
              for \eqn{\delta/\sigma_D} }
  \item{eps2}{right-hand limit of the hypothetical equivalence range for \eqn{\delta/\sigma_D} }
  \item{tol}{tolerance for the error induced through truncating the range of integration on the right}
  \item{PPOST}{posterior probability of the set of all \eqn{(\delta,\sigma_D)} such that
                  \eqn{-\varepsilon_1 < \delta/\sigma_D < \varepsilon_2} }
 }
\references{
Wellek S: Testing statistical hypotheses of equivalence and 
noninferiority. Second edition. Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 3.2.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
postmys(23,0.16,3.99,0.5,0.5,1e-6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{paired data}
\keyword{normal distribution}
\keyword{standardized difference of means}
\keyword{Bayesian test}
\keyword{Jeffreys prior}
\keyword{Gauss-Legendre quadrature}
