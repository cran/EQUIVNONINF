\name{sgnrk}
\alias{sgnrk}
\title{Signed rank test for equivalence of an arbitrary continuous distribution of the intraindividual 
differences in terms of the probability of a positive sign of
a Walsh average: test statistic and critical upper bound}
\description{Implementation of the paired-data analogue of the Mann-Whitney-Wilcoxon test for
equivalence of continuous distributions. The continuity assumption relates to the 
intraindividual differences \eqn{D_i}. For details see Wellek S (2010) Testing statistical                     
hypotheses of equivalence and noninferiority. Second edition,\eqn{\S}{Par.} 5.4.
}
\usage{
sgnrk(alpha,n,qpl1,qpl2,d) 
}
\arguments{
   \item{alpha}{significance level}
   \item{n}{sample size}
   \item{qpl1}{lower equivalence limit \eqn{q_+^{\prime}} to the target functional \eqn{q_+}}
   \item{qpl2}{upper equivalence limit \eqn{q_+^{\prime\prime}} to the target functional \eqn{q_+}}
   \item{d}{row vector with the intraindividual differences for all \eqn{n} pairs as components}
 }
\details{\eqn{q_+} is the probability of getting a positive sign of the so-called Walsh-average 
         of a pair of within-subject differences and can be viewed as a natural paired-observations
         analogue of the Mann-Whitney functional \eqn{\pi_+ = P[X>Y]}.}
\value{
   \item{alpha}{significance level}
   \item{n}{sample size}
   \item{qpl1}{lower equivalence limit \eqn{q_+^{\prime}} to the target functional \eqn{q_+}}
   \item{qpl2}{upper equivalence limit \eqn{q_+^{\prime\prime}} to the target functional \eqn{q_+}}
  \item{U_pl}{observed value of the \eqn{U}-statistics estimator of \eqn{q_+}}
  \item{SIGMAH}{square root of the estimated asymtotic variance of \eqn{U_+}}
  \item{CRIT}{upper critical bound to \eqn{\big|U_+-\big(q_+^{\prime}+q_+^{\prime\prime}\big)/2\big|/\hat{\sigma}}}
  \item{REJ}{indicator of a positive [=1] vs negative [=0] rejection decision to be taken with
             the data under analysis}
  }
\references{
Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 5.4.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
d <- c(-0.5,0.333,0.667,1.333,1.5,-2.0,-1.0,-0.167,1.667,0.833,-2.167,-1.833,
       4.5,-7.5,2.667,3.333,-4.167,5.667,2.333,-2.5)
sgnrk(0.05,20,0.2398,0.7602,d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{nonparametric paired-sample equivalence test}
\keyword{continuous observations}
\keyword{Walsh averages}
\keyword{asymptotic normality}


