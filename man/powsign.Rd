\name{powsign}
\alias{powsign}
\title{Nonconditional power of the UMPU sign test for equivalence and its 
       nonrandomized counterpart}
\description{The program computes for each possible value of the number \eqn{n_0} of 
             zero observations the power conditional on \eqn{N_0 = n_0} and averages
             these conditional power values with respect to the distribution of \eqn{N_0}.
             Equivalence is defined in terms of the logarithm of the ratio \eqn{p_+/p_-}, where
             \eqn{p_+} and \eqn{p_-} denotes the probability of obtaining a positive and negative
             sign, respectively.}
\usage{
 powsign(alpha,n,eps1,eps2,poa)
}
\arguments{
  \item{alpha}{significance level}
  \item{n}{sample size}
  \item{eps1}{absolute value of the lower limit of the hypothetical equivalence range for
              \eqn{\log(p_+/p_-)}.}
  \item{eps2}{upper limit of the hypothetical equivalence range for \eqn{\log(p_+/p_-)}.}
  \item{poa}{probability of a tie under the alternative of interest}
}
\value{
  \item{alpha}{significance level}
  \item{n}{sample size}
  \item{eps1}{absolute value of the lower limit of the hypothetical equivalence range for
              \eqn{\log(p_+/p_-)}.}
  \item{eps2}{upper limit of the hypothetical equivalence range for \eqn{\log(p_+/p_-)}.}
  \item{poa}{probability of a tie under the alternative of interest}
   \item{POWNONRD}{power of the nonrandomized version of the test against the alternative
                   \eqn{p_+ = p_- = (1-p_\circ)/2}}
   \item{POW}{power of the randomized UMPU test against the alternative
                   \eqn{p_+ = p_- = (1-p_\circ)/2}} 
 }
\references{
Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 5.1.
}
\author{
Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>
}
\note{A special case of the test whose power is computed by this program, is the exact conditional 
equivalence test for the McNemar setting (cf. Wellek 2010, pp. 76-77).
}
\examples{    
powsign(0.06580,50,0.847298,0.847298,0.26)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{paired samples}
\keyword{binary observations}
\keyword{McNemar setting}
\keyword{equivalence}
\keyword{exact conditional test}
\keyword{UMPU test}


