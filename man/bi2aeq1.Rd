\name{bi2aeq1}
\alias{bi2aeq1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Power of the exact Fisher type test for equivalence}
\description{The function computes exact values of the power of the randomized UMPU test for 
equivalence in the strict (i.e. two-sided) sense of two binomial distributions and the 
conservative nonrandomized version of that test. It is assumed that the samples 
being available from both distributions are independent.
}
\usage{
bi2aeq1(m,n,rho1,rho2,alpha,p1,p2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{size of Sample 1}
  \item{n}{size of Sample 2}
  \item{rho1}{lower limit of the hypothetical equivalence range for the odds ratio}
  \item{rho2}{upper limit of the hypothetical equivalence range for the odds ratio}
  \item{alpha}{significance level}
  \item{p1}{true success rate in Population 1}
  \item{p2}{true success rate in Population 2}
}
\value{
  \item{m}{size of Sample 1}
  \item{n}{size of Sample 2}
  \item{rho1}{lower limit of the hypothetical equivalence range for the odds ratio}
  \item{rho2}{upper limit of the hypothetical equivalence range for the odds ratio}
  \item{alpha}{significance level}
  \item{p1}{true success rate in Population 1}
  \item{p2}{true success rate in Population 2}
  \item{POWNR}{Power of the nonrandomized version of the test}
  \item{POW}{Power of the randomized UMPU test}
}
\references{Wellek S: Testing statistical hypotheses of equivalence and noninferiority. 
        Second edition. Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 6.6.4.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
bi2aeq1(302,302,0.6667,1.5,0.05,0.5,0.5)
  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{binomial two-sample problem}
\keyword{equivalence}
\keyword{odds ratio}
\keyword{exact Fisher-type test}
\keyword{power}
