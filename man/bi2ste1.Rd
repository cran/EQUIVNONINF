\name{bi2ste1}
\alias{bi2ste1}
\title{Power of the exact Fisher type test for noninferiority}
\description{The function computes exact values of the power of the randomized UMPU test for one-sided 
             equivalence of two binomial distributions and its conservative nonrandomized version. 
             It is assumed that the samples being available from both distributions are independent.
}
\usage{
bi2ste1(m, n, eps, alpha, p1, p2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{size of Sample 1}
  \item{n}{size of Sample 2}
  \item{eps}{noninferiority margin to the odds ratio \eqn{\varrho}, defined to be the maximum 
             acceptable deviation of the true value of \eqn{\varrho} from unity}
  \item{alpha}{significance level}
  \item{p1}{success rate in Population 1}
  \item{p2}{success rate in Population 2}
}
\value{
  \item{m}{size of Sample 1}
  \item{n}{size of Sample 2}
  \item{eps}{noninferiority margin to the odds ratio \eqn{\varrho}, defined to be the maximum 
             acceptable deviation of the true value of \eqn{\varrho} from unity}
  \item{alpha}{significance level}
  \item{p1}{success rate in Population 1}
  \item{p2}{success rate in Population 2}
\item{POWNR}{power of the nonrandomized version of the test}
\item{POW}{power of the randomized UMPU test}
%% ...
}
\references{Wellek S: Testing statistical hypotheses of equivalence and noninferiority. 
            Second edition. Boca Raton: Chapman & Hall/CRC Press, 2010, 6.6.1.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
bi2ste1(106,107,0.5,0.05,0.9245,0.9065)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{binomial two-sample problem}
\keyword{noninferiority}
\keyword{odds ratio}
\keyword{exact Fisher-type test}
\keyword{power}
