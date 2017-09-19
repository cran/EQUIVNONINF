\name{bi2by_ni_OR}
\alias{bi2by_ni_OR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Objective Bayesian test for noninferiority in the two-sample setting with
binary data and the odds ratio as the parameter of interest
}
\description{Implementation of the construction described on pp. 179--181 of
Wellek S (2010) Testing statistical hypotheses of equivalence and 
noninferiority. Second edition.
}
\usage{
bi2by_ni_OR(N1,N2,EPS,SW,NSUB,ALPHA,MAXH)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N1}{size of sample 1}
  \item{N2}{size of sample 2}
  \item{EPS}{noninferiority margin to the deviation of the odds ratio from unity}
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
  \item{N1}{size of sample 1}
  \item{N2}{size of sample 2}
  \item{EPS}{noninferiority margin to the deviation of the odds ratio from unity}
  \item{NSUB}{number of subintervals for partitioning the range of integration}
  \item{SW}{width of the search grid for determining the maximum of the rejection probability
           on the common boundary of the hypotheses} 
  \item{ALPHA0}{result of the search for the largest admissible nominal level}
  \item{SIZE0}{size of the critical region corresponding to \eqn{\alpha_0}{alpha_0}}
  \item{SIZE_UNC}{size of the critical region of the test at uncorrected nominal level}
}
\references{Wellek S: Statistical methods for the analysis of two-arm non-inferiority trials
            with binary outcomes. Biometrical Journal 47 (2005), 48--61.

        Wellek S: Testing statistical hypotheses of equivalence and 
        noninferiority. Second edition. Boca Raton: 
        Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 6.6.2.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
bi2by_ni_OR(10,10,1/3,.0005,10,.05,12)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{parallel-group design with binary data}
\keyword{noninferiority}
\keyword{odds ratio}
\keyword{Bayesian test}
\keyword{product of Jeffrey's priors}
\keyword{Gauss-Legendre quadrature}