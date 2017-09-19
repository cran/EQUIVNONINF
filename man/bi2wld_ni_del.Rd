\name{bi2wld_ni_del}
\alias{bi2wld_ni_del}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute corrected nominal levels for the Wald type
(asymptotic) test for one-sided equivalence of two binomial distributions
with respect to the difference of success rates
}
\description{Implementation of the construction described on pp. 183-5 of
Wellek S (2010) Testing statistical hypotheses of equivalence and 
                noninferiority. Second edition.
}
\usage{
bi2wld_ni_del(N1,N2,EPS,SW,ALPHA,MAXH)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N1}{size of Sample 1}
  \item{N2}{size of Sample 2}
  \item{EPS}{noninferiority margin to the difference of success probabilities}
  \item{SW}{width of the search grid for determining the maximum of the rejection probability
           on the common boundary of the hypotheses} 
  \item{ALPHA}{target significance level}
  \item{MAXH}{maximum number of interval-halving steps}
}
\details{The program computes the largest nominal significance level
to be used for determining the critical lower bound to the Wald-type statistic for the
problem of testing \eqn{H:p_1 \le p_2 - \varepsilon} versus \eqn{K: p_1 < p_2 - \varepsilon}.}
\value{
  \item{N1}{size of Sample 1}
  \item{N2}{size of Sample 2}
  \item{EPS}{noninferiority margin to the difference of success probabilities}
  \item{SW}{width of the search grid for determining the maximum of the rejection probability
           on the common boundary of the hypotheses} 
  \item{ALPHA}{target significance level}
  \item{MAXH}{maximum number of interval-halving steps}
\item{ALPHA0}{corrected nominal level}
\item{SIZE0}{size of the critical region of the test at nominal level ALPHA0}
\item{SIZE_UNC}{size of the test at uncorrected nominal level ALPHA}
\item{ERR_IND}{indicator taking value 1 when it occurs that the sufficient condition
allowing one to restrict the search for the maximum of the rejection probability 
under the null hypothesis to its boundary, fails to be satisfied; otherwise the
indicator retains its default value 0.}
}
\references{Wellek S: Testing statistical hypotheses of equivalence and 
        noninferiority. Second edition. Boca Raton: 
        Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 6.6.3.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
bi2wld_ni_del(25,25,.10,.01,.05,10) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Wald-type}
\keyword{nominal level}
