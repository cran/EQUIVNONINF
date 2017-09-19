\name{po_pbibe}
\alias{po_pbibe}
\title{Bayesian posterior probability of the alternative
hypothesis of probability-based individual bioequivalence (PBIBE)}
\description{Implementation of the algorithm presented in \eqn{\S}{Par.} 10.3.3 of
Wellek S (2010) Testing statistical hypotheses of equivalence and 
noninferiority. Second edition.}
\usage{
po_pbibe(n,eps,pio,zq,s,tol,sw,ihmax)
}
\arguments{
  \item{n}{sample size}
  \item{eps}{equivalence margin to an individual log-bioavailability ratio}
  \item{pio}{prespecified lower bound to the probability of obtaining an individual
            log-bioavailability ratio falling in the equivalence range \eqn{(-\varepsilon,\varepsilon)}}
  \item{zq}{mean log-bioavailability ratio observed in the sample under analysis} 
  \item{s}{square root of the sample variance of the log-bioavailability ratios} 
  \item{tol}{maximum numerical error allowed for transforming the hypothesis of PBIBE 
             into a region in the parameter space of the log-normal distribution assumed 
             to underlie the given sample of individual bioavailability ratios}
  \item{sw}{step width used in the numerical procedure yielding results at a level of accuracy
            specified by the value chosen for tol}
  \item{ihmax}{maximum number of interval halving steps to be carried out in finding the
               region specified in the parameter space according to the criterion  of PBIBE}
}
\details{The program uses 96-point Gauss-Legendre quadrature.}
\value{
  \item{n}{sample size}
  \item{eps}{equivalence margin to an individual log-bioavailability ratio}
  \item{pio}{prespecified lower bound to the probability of obtaining an individual
            log-bioavailability ratio falling in the equivalence range \eqn{(-\varepsilon,\varepsilon)}}
  \item{zq}{mean log-bioavailability ratio observed in the sample under analysis} 
  \item{s}{square root of the sample variance of the log-bioavailability ratios} 
  \item{tol}{maximum numerical error allowed for transforming the hypothesis of PBIBE 
             into a region in the parameter space of the log-normal distribution assumed 
             to underlie the given sample of individual bioavailability ratios}
  \item{sw}{step width used in the numerical procedure yielding results at a level of accuracy
            specified by the value chosen for tol}
  \item{ihmax}{maximum number of interval halving steps to be carried out in finding the
               region specified in the parameter space according to the criterion  of PBIBE}
  \item{PO_PBIBE}{posterior probability of the alternative hypothesis of PBIBE}
 }
\references{
Wellek S: Bayesian construction of an improved parametric test for probability-based 
individual bioequivalence. Biometrical Journal 42 (2000), 1039-52.

Wellek S: Testing statistical hypotheses of equivalence and 
noninferiority. Second edition. Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 10.3.3.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
po_pbibe(20,0.25,0.75,0.17451,0.04169, 10e-10,0.01,100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{individual bioequivalence}
\keyword{lognormal distribution}
\keyword{Bayesian test}
\keyword{Gauss-Legendre quadrature}
