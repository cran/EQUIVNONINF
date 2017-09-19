\name{mcnby_ni}
\alias{mcnby_ni}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian test for noninferiority in the McNemar setting with the
difference of proportions as the parameter of interest}
\description{
The program determines through iteration the largest nominal
level \eqn{\alpha_0} such that comparing the posterior probability
of the alternative hypothesis \eqn{K_1: \delta > -\delta_0} to the lower
bound \eqn{1-\alpha_0} generates a critical region whose size does not exceed
the target significance level \eqn{\alpha}. In addition, exact values of the
power against specific parameter configurations with \eqn{\delta = 0} are output.
}
\usage{
mcnby_ni(N,DEL0,K1,K2,K3,NSUB,SW,ALPHA,MAXH)
}
\arguments{
  \item{N}{sample size}
  \item{DEL0}{noninferiority margin to the difference of the parameters of the marginal
              binomial distributions under comparison}
  \item{K1}{Parameter 1 of the Dirichlet prior for the family of trinomial distributions}
  \item{K2}{Parameter 2 of the Dirichlet prior for the family of trinomial distributions}
  \item{K3}{Parameter 3 of the Dirichlet prior for the family of trinomial distributions}
  \item{NSUB}{number of subintervals for partitioning the range of integration}
  \item{SW}{width of the search grid for determining the maximum of the rejection 
            probability on the common boundary of the hypotheses}
  \item{ALPHA}{target significance level}
  \item{MAXH}{maximum number of interval halving steps to be carried out in finding the
              maximally raised nominal level}
  
}
\details{The program uses 96-point Gauss-Legendre quadrature on each of
the NSUB intervals into which the range of integration is partitioned.}
\value{
  \item{N}{sample size}
  \item{DEL0}{noninferiority margin to the difference of the parameters of the marginal
              binomial distributions under comparison}
  \item{K1}{Parameter 1 of the Dirichlet prior for the family of trinomial distributions}
  \item{K2}{Parameter 2 of the Dirichlet prior for the family of trinomial distributions}
  \item{K3}{Parameter 3 of the Dirichlet prior for the family of trinomial distributions}
  \item{NSUB}{number of subintervals for partitioning the range of integration}
  \item{SW}{width of the search grid for determining the maximum of the rejection 
            probability on the common boundary of the hypotheses}
  \item{ALPHA}{target significance level}
  \item{MAXH}{maximum number of interval halving steps to be carried out in finding the
              maximally raised nominal level}
  \item{ALPHA0}{result of the search for the largest admissible nominal level}
  \item{SIZE0}{size of the critical region corresponding to \eqn{\alpha_0}}
  \item{SIZE_UNC}{size of the critical region of test at uncorrected nominal level \eqn{\alpha}}
  \item{POW}{power against 7 different parameter configurations with \eqn{\delta =0}}
}
\references{Wellek S: Testing statistical hypotheses of equivalence and 
        noninferiority. Second edition. Boca Raton: 
        Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 5.2.3.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
mcnby_ni(25,.10,.5,.5,.5,10,.05,.05,5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{matched-pair design with binary data}
\keyword{noninferiority}
\keyword{difference of proportions}
\keyword{Bayesian test}
\keyword{Dirichlet prior}
\keyword{Gauss-Legendre quadrature}