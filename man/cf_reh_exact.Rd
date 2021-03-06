\name{cf_reh_exact}
\alias{cf_reh_exact}
\title{Exact confidence bounds to the relative excess heterozygosity (REH) exhibited
by a SNP genotype distribution
}
\description{Implementation of the interval estimation procedure described on 
   pp. 305-6 of Wellek S (2010) Testing statistical hypotheses of equivalence and 
   noninferiority. Second edition.
}
\usage{
cf_reh_exact(X1,X2,X3,alpha,SW,TOL,ITMAX)
}
\arguments{
  \item{X1}{count of homozygotes of the first kind [\eqn{\leftrightarrow}{<->} genotype AA]}
  \item{X2}{count of heterozygotes [\eqn{\leftrightarrow}{<->} genotype AB]}
  \item{X3}{count of homozygotes of the second kind [\eqn{\leftrightarrow}{<->} genotype BB]}
  \item{alpha}{1 - confidence level} 
  \item{SW}{width of the search grid for determining an interval covering the 
            parameter point at which the conditional distribution function takes
            value \eqn{\alpha} and \eqn{1-\alpha}, respectively} 
  \item{TOL}{numerical tolerance to the deviation between the computed confidence limits
             and their exact values}
  \item{ITMAX}{maximum number of interval-halving steps}
}
\details{The program exploits the structure of the family of all genotype distributions,
         which is 2-parameter exponential with \eqn{\log(REH)} as one of these parameters.}
\value{
  \item{X1}{count of homozygotes of the first kind [\eqn{\leftrightarrow}{<->} genotype AA]}
  \item{X2}{count of heterozygotes [\eqn{\leftrightarrow}{<->} genotype AB]}
  \item{X3}{count of homozygotes of the second kind [\eqn{\leftrightarrow}{<->} genotype BB]}
  \item{alpha}{1 - confidence level} 
  \item{SW}{width of the search grid for determining an interval covering the 
            parameter point at which the conditional distribution function takes
            value \eqn{\alpha} and \eqn{1-\alpha}, respectively} 
  \item{TOL}{numerical tolerance to the deviation between the computed confidence limits
             and their exact values}
  \item{ITMAX}{maximum number of interval-halving steps}
\item{C_l_exact}{exact conditional lower \eqn{(1-\alpha)}-confidence bound to REH}
\item{C_r_exact}{exact conditional upper \eqn{(1-\alpha)}-confidence bound to REH}
}
\references{
Wellek S, Goddard KAB, Ziegler A: A confidence-limit-based approach
to the assessment of Hardy-Weinberg equilibrium. 
Biometrical Journal 52 (2010), 253-270.

Wellek S: Testing statistical hypotheses of equivalence and 
noninferiority. Second edition. Boca Raton: 
Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 9.4.3.}

\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{               
cf_reh_exact(34,118,96,.05,.1,1E-4,25)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{genetic association study}
\keyword{Hardy-Weinberg equilibrium}
\keyword{goodness of fit}
\keyword{equivalence test}
\keyword{confidence interval inclusion}
\keyword{relative excess heterozygosity}
\keyword{exact conditional confidence limit}


