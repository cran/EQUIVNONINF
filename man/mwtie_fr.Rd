\name{mwtie_fr}
\alias{mwtie_fr}
\title{Analogue of mwtie_xy for settings with grouped data}
\description{Implementation of the asymptotically distribution-free test for 
equivalence of discrete distributions from which grouped data are obtained.
Hypothesis formulation is in terms of the Mann-Whitney-Wilcoxon functional
generalized to the case that ties between observations from different distributions may
occur with positive probability. For details see Wellek S (2010) Testing statistical hypotheses 
of equivalence and noninferiority. Second edition, p.155.}
\usage{
mwtie_fr(k,alpha,m,n,eps1_,eps2_,x,y) 
}
\arguments{
   \item{k}{total number of grouped values which can be distinguished in the pooled sample}
   \item{alpha}{significance level}
   \item{m}{size of Sample 1} 
   \item{n}{size of Sample 2}
   \item{eps1_}{absolute value of the left-hand limit of the hypothetical equivalence range for
                \eqn{\pi_+/(1-\pi_0) - 1/2}}
   \item{eps2_}{right-hand limit of the hypothetical equivalence range for \eqn{\pi_+/(1-\pi_0) - 1/2}}
   \item{x}{row vector with the \eqn{m} observations making up Sample1 as components}
   \item{y}{row vector with the \eqn{n} observations making up Sample2 as components}
 }
\details{Notation: \eqn{\pi_+} and \eqn{\pi_0} stands for the functional defined by \eqn{\pi_+ = P[X>Y]} and
\eqn{\pi_0 = P[X=Y]}, respectively, 
with \eqn{X\sim F \equiv} cdf of Population 1 being independent of \eqn{Y\sim G \equiv} cdf of Population 2.
}
\value{
   \item{alpha}{significance level}
   \item{m}{size of Sample 1} 
   \item{n}{size of Sample 2}
   \item{eps1_}{absolute value of the left-hand limit of the hypothetical equivalence range for
                \eqn{\pi_+/(1-\pi_0) - 1/2}}
   \item{eps2_}{right-hand limit of the hypothetical equivalence range for \eqn{\pi_+/(1-\pi_0) - 1/2}}
  \item{WXY_TIE}{observed value of the \eqn{U}-statistics -- based estimator of \eqn{\pi_+/(1-\pi_0)}}
  \item{SIGMAH}{square root of the estimated asymtotic variance of \eqn{W_+/(1-W_0)}}
  \item{CRIT}{upper critical bound to \eqn{|W_+/(1-W_0) - 1/2 - 
               (\varepsilon^\prime_2-\varepsilon^\prime_1)/2|/\hat{\sigma}}}
  \item{REJ}{indicator of a positive [=1] vs negative [=0] rejection decision to be taken with
             the data under analysis}
  }
\references{
Wellek S, Hampel B: A distribution-free two-sample equivalence test allowing for tied 
observations. Biometrical Journal 41 (1999), 171-186.

Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 6.4.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
x <- c(1,1,3,2,2,3,1,1,1,2,1,2,2,2,1,2,1,3,2,1,2,1,1,1,1,1,1,1,1,1,1,1,2,1,3,1,3,2,1,1,
       2,1,2,1,1,2,2,1,2,1,1,1,1,1,2,2,1,2,2,1,3,1,2,1,1,2,2,1,2,2,1,1,1,3,2,1,1,1,2,1,
       3,3,3,1,2,1,2,2,1,1,1,2,2,1,1,2,1,1,2,3,1,3,2,1,1,1,1,2,2,2,1,1,2,2,3,2,1,2,1,1,
       2,2,1,2,2,2,1,1,2,3,2,1,3,2,1,1,1,2,2,2,2,1,2,2,1,1,1,1,2,1,1,1,2,1,2,2,1,2,2,2,
       2,1,1,2,1,2,2,1,1,1,1,3,1,1,2,2,1,1,1,2,2,2,1,2,3,2,2,1,2,1,2,1,1,2,1,2,2,1,1,1,
       2,2,2,2)
y <- c(2,1,2,2,1,1,2,2,2,1,1,2,1,3,3,1,1,1,1,1,1,2,2,3,1,1,1,3,1,1,1,1,1,1,1,2,2,3,2,1,
       2,2,2,1,2,1,1,2,2,1,2,1,1,1,1,2,1,2,1,1,3,1,1,1,2,2,2,1,1,1,1,2,1,2,1,1,2,2,2,2,
       2,1,1,1,3,2,2,2,1,2,3,1,2,1,1,1,2,1,3,3,1,2,2,2,2,2,2,1,2,1,1,1,1,2,2,1,1,1,1,2,
       1,3,1,1,2,1,2,1,2,2,2,1,2,2,2,1,1,1,2,1,2,1,2,1,1,1,2,1,2,2,1,1,1,1,2,2,3,1,3,1,
       1,2,2,2,1,1,1,1,2,1,1,3,2,2,3,1,2,2,1,1,2,1,1,2,1,2,2,1,2,1,2,2,2,1,1,1,1,1,1,1,
       1,1,1,2,1,3,2,2,1,1,1,2,2,1,1,2,1,2,1,2,2,2,1,2,3,1,1,2,1,2,2,1,1,1,1,2,2,2,1,1,
       3,2,1,2,2,2,1,1,1,2,1,2,2,1,2,1,1,2)
mwtie_fr(3,0.05,204,258,0.10,0.10,x,y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{Mann-Whitney functional}
\keyword{nonparametric two-sample equivalence test}
\keyword{tied observations}
\keyword{U-statistics estimators}
\keyword{asymptotic normality}


