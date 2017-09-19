\name{pow_abe}
\alias{pow_abe}
\title{Confidence innterval inclusion test for average bioequivalence:
exact power against an arbitrary specific alternative}
\description{
Evaluation of the integral on the right-hand side of equation (10.11)
of p. 317 of Wellek S (2010) Testing statistical hypotheses of equivalence and 
noninferiority. Second edition}
\usage{
pow_abe(m,n,alpha,del_0,del,sig)
}
\arguments{
  \item{m}{sample size in sequence group T(est)/R(eference)}
  \item{n}{sample size in sequence group R(eference)/T(est)}
  \item{alpha}{significance level} 
  \item{del_0}{equivalence margin to the absolute value of the log-ratio \eqn{\mu^*_T} and \eqn{\mu^*_R} 
              of the formulation effects}
  \item{del}{assumed true value of \eqn{| \log(\mu^*_T/\mu^*_R) |}, with \eqn{0\le\delta < \delta_0} }
  \item{sig}{theoretical standard deviation of the log within-subject bioavailability ratios
             in each sequence group} 
}
\details{The program uses 96-point Gauss-Legendre quadrature.}
\value{
  \item{m}{sample size in sequence group T(est)/R(eference)}
  \item{n}{sample size in sequence group R(eference)/T(est)}
  \item{alpha}{significance level} 
  \item{del_0}{equivalence margin to the absolute value of the log-ratio \eqn{\mu^*_T} and \eqn{\mu^*_R} 
              of the formulation effects}
  \item{del}{assumed true value of \eqn{| \log(\mu^*_T/\mu^*_R) |}, with \eqn{0\le\delta < \delta_0} }
  \item{POW_ABE}{power of the interval inclusion test for average bioequivalence against the
                 specific alternative given by \eqn{(\delta,\sigma)}}
}
\references{
Wellek S: Testing statistical hypotheses of equivalence and 
noninferiority. Second edition. Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 10.2.1.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
pow_abe(12,13,0.05,log(1.25),log(1.25)/2,0.175624)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{average bioequivalence}
\keyword{confidence interval inclusion test}
\keyword{two-period crossover}
\keyword{lognormal distribution}
