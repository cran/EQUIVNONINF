\name{bi2rlv2}
\alias{bi2rlv2}
\title{Sample sizes for the exact Fisher type test for relevant differences}
\description{The function computes minimum sample sizes required in the randomized UMPU test for 
relevant differences between two binomial distributions with respect to the odds ratio. Computation is done 
under the side condition  that the ratio \eqn{m/n} has some predefined value \eqn{\lambda}.
}
\usage{
bi2rlv2(rho1,rho2,alpha,p1,p2,beta,qlambd)
}
\arguments{
  \item{rho1}{lower limit of the hypothetical equivalence range for the odds ratio}
  \item{rho2}{upper limit of the hypothetical equivalence range for the odds ratio}
  \item{alpha}{significance level}
  \item{p1}{true success rate in Population 1}
  \item{p2}{true success rate in Population 2}
  \item{beta}{target value of power}
  \item{qlambd}{sample size ratio \eqn{m/n}}
}
\value{
  \item{rho1}{lower limit of the hypothetical equivalence range for the odds ratio}
  \item{rho2}{upper limit of the hypothetical equivalence range for the odds ratio}
  \item{alpha}{significance level}
  \item{p1}{true success rate in Population 1}
  \item{p2}{true success rate in Population 2}
  \item{beta}{target value of power}
  \item{qlambd}{sample size ratio \eqn{m/n}}
  \item{M}{minimum size of Sample 1}
  \item{N}{minimum size of Sample 2}
  \item{POW}{power of the randomized UMPU test attained with the computed values of m, n}
}
\references{Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 11.3.3.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{                                     
bi2rlv2(.6667,1.5,.05,.70,.50,.80,2.0)  }  
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{binomial two-sample problem}
\keyword{relevant differences}
\keyword{odds ratio}
\keyword{exact Fisher-type test}
\keyword{sample size}