\name{bi2ste3}
\alias{bi2ste3}
\title{Determination of a maximally raised nominal significance level
for the nonrandomized version of the exact Fisher type test for noninferiority
}
\description{The objective is to raise the nominal significance level as far as
possible without exceeding the target significance level in the nonrandomized
version of the test. The approach goes back to R.D. Boschloo (1970) who used the
same technique for reducing the conservatism of the traditional nonrandomized 
Fisher test for superiority. }
\usage{
bi2ste3(m, n, eps, alpha, sw, tolrd, tol, maxh)
}
\arguments{
  \item{m}{size of Sample 1}
  \item{n}{size of Sample 2}
  \item{eps}{noninferiority margin to the odds ratio \eqn{\varrho}, defined to be the maximum     
             acceptable deviation of the true value of \eqn{\varrho} from unity}
 \item{alpha}{target significance level}
 \item{sw}{width of the search grid for determining the maximum of the rejection probability
           on the common boundary of the hypotheses}
  \item{tolrd}{horizontal distance from 0 and 1, respectively, of the left- and right-most boundary point 
               to be included in the search grid}
  \item{tol}{upper bound to the absolute difference between size and target level below
             which the search for a corrected nominal level terminates}
  \item{maxh}{maximum number of interval-halving steps to be carried out in finding the
              maximally raised nominal level}
}
\details{It should be noted that, as the function of the nominal level, the size of the
        nonrandomized test is piecewise constant. Accordingly, there is a nondegenerate 
        interval of "candidate" nominal levels serving the purpose. The limits of such
        an interval can be read from the output.}

\value{
   \item{m}{size of Sample 1}
   \item{n}{size of Sample 2}
   \item{eps}{noninferiority margin to the odds ratio \eqn{\varrho}, defined to be the maximum     
             acceptable deviation of the true value of \eqn{\varrho} from unity}
   \item{alpha}{target significance level}
   \item{sw}{width of the search grid for determining the maximum of the rejection probability
           on the common boundary of the hypotheses}
   \item{tolrd}{horizontal distance from 0 and 1, respectively, of the left- and right-most boundary point 
               to be included in the search grid}
   \item{tol}{upper bound to the absolute difference between size and target level below
             which the search for a corrected nominal level terminates}
   \item{maxh}{maximum number of interval-halving steps to be carried out in finding the
              maximally raised nominal level}
   \item{ALPH_0}{current trial value of the raised nominal level searched for}
   \item{NHST}{number of interval-halving steps performed up to now}
   \item{SIZE}{size of the critical region corresponding to \eqn{\alpha_0}{alpha_0}}
}
\references{
Boschloo RD: Raised conditional level of significance for the 2 x 2-
table when testing the equality of two probabilities. Statistica Neerlandica 24 (1970), 
1-35.

Wellek S: Testing statistical hypotheses of equivalence and noninferiority. 
          Second edition. Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}6.6.2.
}
\author{
Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>
}
\examples{
 bi2ste3(50, 50, 1/3, 0.05, 0.05, 1e-10, 1e-8, 10) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{binomial two-sample problem}
\keyword{noninferiority}
\keyword{odds ratio}
\keyword{exact Fisher-type test nonrandomized version}
\keyword{maximally raised nominal level}

