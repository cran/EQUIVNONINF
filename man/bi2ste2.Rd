\name{bi2ste2}
\alias{bi2ste2}
\title{Sample sizes for the exact Fisher type test for noninferiority}
\description{Sample sizes for the exact Fisher type test for noninferiority
}
\usage{
bi2ste2(eps, alpha, p1, p2, bet, qlambd)}
\arguments{
  \item{eps}{noninferiority margin to the odds ratio}
  \item{alpha}{significance level}
  \item{p1}{success rate in Population 1}
  \item{p2}{success rate in Population 2}
  \item{bet}{target power value}
  \item{qlambd}{sample size ratio \eqn{m/n} }
}
\details{The program computes the smallest sample sizes \eqn{m},\eqn{n} satisfying
\eqn{m/n = \lambda} required for ensuring that the power of the randomized UMPU test does not 
  fall below \eqn{\beta}. }
\value{
  \item{eps}{noninferiority margin to the odds ratio}
  \item{alpha}{significance level}
  \item{p1}{success rate in Population 1}
  \item{p2}{success rate in Population 2}
  \item{bet}{target power value}
  \item{qlambd}{sample size ratio \eqn{m/n} }
\item{M}{minimum size of Sample 1}
\item{N}{minimum size of Sample 2}
\item{POW}{power of the randomized UMPU test attained with the computed values of m, n}
}
\references{Wellek S: Testing statistical hypotheses of equivalence and 
        noninferiority. Second edition. Boca Raton: 
        Chapman & Hall/CRC Press, 2010, 6.6.1.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
bi2ste2(0.5,0.05,0.9245,0.9065,0.80,1.00)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{binomial two-sample problem}
\keyword{noninferiority}
\keyword{odds ratio}
\keyword{exact Fisher-type test}
\keyword{sample size}
