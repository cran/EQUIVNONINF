\name{mawi}
\alias{mawi}
\title{Mann-Whitney test for equivalence of two continuous distributions of arbitrary shape:
test statistic and critical upper bound}
\description{Implementation of the asymptotically distribution-free test for 
equivalence of two continuous distributions in terms of the Mann-Whitney-Wilcoxon functional. 
For details see Wellek S (2010) Testing statistical hypotheses of equivalence and 
noninferiority. Second edition, \eqn{\S}{Par.} 6.2.
}
\usage{
mawi(alpha,m,n,eps1_,eps2_,x,y) 
}
\arguments{
   \item{alpha}{significance level}
   \item{m}{size of Sample 1} 
   \item{n}{size of Sample 2}
   \item{eps1_}{absolute value of the left-hand limit of the hypothetical equivalence range for
                \eqn{\pi_+ - 1/2}}
   \item{eps2_}{right-hand limit of the hypothetical equivalence range for \eqn{\pi_+ - 1/2}}
   \item{x}{row vector with the \eqn{m} observations making up Sample1 as components}
   \item{y}{row vector with the \eqn{n} observations making up Sample2 as components}
 }
\details{Notation: \eqn{\pi_+} stands for the Mann-Whitney functional defined by \eqn{\pi_+ = P[X>Y]},
with \eqn{X\sim F \equiv} cdf of Population 1 being independent of \eqn{Y\sim G \equiv} cdf of Population 2.
}
\value{
   \item{alpha}{significance level}
   \item{m}{size of Sample 1} 
   \item{n}{size of Sample 2}
   \item{eps1_}{absolute value of the left-hand limit of the hypothetical equivalence range for
                \eqn{\pi_+ - 1/2}}
   \item{eps2_}{right-hand limit of the hypothetical equivalence range for \eqn{\pi_+ - 1/2}}
  \item{W+}{observed value of the \eqn{U}-statistics estimator of \eqn{\pi_+}}
  \item{SIGMAH}{square root of the estimated asymtotic variance of \eqn{W_+}}
  \item{CRIT}{upper critical bound to \eqn{|W_+ - 1/2 - 
               (\varepsilon^\prime_2-\varepsilon^\prime_1)/2|/\hat{\sigma}}}
  \item{REJ}{indicator of a positive [=1] vs negative [=0] rejection decision to be taken with
             the data under analysis}
  }
\references{
Wellek S: A new approach to equivalence assessment in standard comparative bioavailability trials
by means of the Mann-Whitney statistic. Biometrical Journal 38 (1996), 695-710.

Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 6.2.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
x <- c(10.3,11.3,2.0,-6.1,6.2,6.8,3.7,-3.3,-3.6,-3.5,13.7,12.6)
y <- c(3.3,17.7,6.7,11.1,-5.8,6.9,5.8,3.0,6.0,3.5,18.7,9.6)
mawi(0.05,12,12,0.1382,0.2602,x,y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{Mann-Whitney functional}
\keyword{nonparametric two-sample equivalence test}
\keyword{continuous observations}
\keyword{asymptotic normality}


