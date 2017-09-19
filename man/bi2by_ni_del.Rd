\name{bi2by_ni_del}
\alias{bi2by_ni_del}
\title{Objective Bayesian test for noninferiority in the two-sample setting with
binary data and the difference of the two proportions as the parameter of interest
}
\description{Implementation of the construction described on pp. 185-6 of
Wellek S (2010) Testing statistical hypotheses of equivalence and 
noninferiority. Second edition.
}
\usage{
bi2by_ni_del(N1,N2,EPS,SW,NSUB,ALPHA,MAXH)
}
\arguments{
  \item{N1}{size of Sample 1}
  \item{N2}{size of sample 2}
  \item{EPS}{noninferiority margin to the difference of success probabilities}
  \item{SW}{width of the search grid for determining the maximum of the rejection probability
           on the common boundary of the hypotheses} 
  \item{NSUB}{number of subintervals for partitioning the range of integration}  
  \item{ALPHA}{target significance level}
  \item{MAXH}{maximum number of interval halving steps to be carried out in finding the
              maximally admissible nominal level}
}
\details{The program uses 96-point Gauss-Legendre quadrature on each of
the  NSUB intervals into which the range of integration is partitioned.}
\value{
  \item{N1}{size of Sample 1}
  \item{N2}{size of sample 2}
  \item{EPS}{noninferiority margin to the difference of success probabilities}
  \item{NSUB}{number of subintervals for partitioning the range of integration}
  \item{SW}{width of the search grid for determining the maximum of the rejection probability
           on the common boundary of the hypotheses}
  \item{ALPHA0}{result of the search for the largest admissible nominal level}
  \item{SIZE0}{size of the critical region corresponding to \eqn{\alpha_0}{alpha_0}}
  \item{SIZE_UNC}{size of the critical region of the test at uncorrected nominal level}
}
\references{Wellek S: Statistical methods for the analysis of two-armed non-inferiority 
            trials with binary outcomes. Biometrical Journal 47 (2005), 48--61.
			
Wellek S: Testing statistical hypotheses of equivalence and 
        noninferiority. Second edition. Boca Raton: 
        Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 6.6.3.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
bi2by_ni_del(20,20,.10,.01,10,.05,5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{parallel-group design with binary data}
\keyword{noninferiority}
\keyword{difference of proportions}
\keyword{Bayesian test}
\keyword{product of Jeffrey's priors}
\keyword{Gauss-Legendre quadrature}
