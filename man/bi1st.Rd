\name{bi1st}
\alias{bi1st}
\title{Critical constants and power  of the UMP test for equivalence of a single binomial proportion 
to some given reference value}
\description{The function computes the critical constants defining the uniformly most powerful (randomized) test 
for the problem \eqn{p \le p_1} or \eqn{p \ge p_2} versus
\eqn{p_1 < p < p_2}, with \eqn{p} denoting the parameter of
a binomial distribution from which a single sample of size
\eqn{n} is available. In the output, one also finds the power 
against the alternative that the true value of \eqn{p} falls on the 
midpoint of the hypothetical equivalence interval \eqn{(p_1 , p_2).} 
}
\usage{
bi1st(alpha,n,P1,P2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{alpha}{significance level}
   \item{n}{sample size}
  \item{P1}{lower limit of the hypothetical equivalence
            range for the binomial parameter \eqn{p}}
  \item{P2}{upper limit of the hypothetical equivalence
            range for \eqn{p}}
}
\value{
  \item{alpha}{significance level}
  \item{n}{sample size}
  \item{P1}{lower limit of the hypothetical equivalence
            range for the binomial parameter \eqn{p}}
  \item{P2}{upper limit of the hypothetical equivalence
            range for \eqn{p}}
\item{C1}{left-hand limit of the critical interval for
          the observed number \eqn{X} of successes}
\item{C2}{right-hand limit of the critical interval for
          \eqn{X}}
\item{GAM1}{probability of rejecting the null hypothesis 
            when it turns out that \eqn{X=C_1}}
\item{GAM2}{probability of rejecting the null hypothesis 
            for \eqn{X=C_2}}

\item{POWNONRD}{Power of the nonrandomized version of the test against the alternative \eqn{p = (p_1+p_2)/2}}
\item{POW}{Power of the randomized UMP test against the        
           alternative \eqn{p = (p_1+p_2)/2}}
}
\references{Wellek S: Testing statistical hypotheses of equivalence and noninferiority. 
            Second edition. Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 4.3.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
        Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
bi1st(.05,273,.65,.75)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{binomial one-sample problem}
\keyword{equivalence}
\keyword{UMP test}
\keyword{randomized test}
\keyword{power}
