\name{mcnemasc}
\alias{mcnemasc}
\title{Determination of a corrected nominal significance level for the asymptotic test for equivalence
of two paired binomial proportions with respect to the difference of their expectations (McNemar
setting)
}
\description{The program computes the largest nominal significance level 
which can be substituted for the target level \eqn{\alpha} without making the exact 
size of the asymptotic testing procedure larger than \eqn{\alpha}.                  }
\usage{
 mcnemasc(alpha,n,del0,sw,tol,maxh)
}
\arguments{
  \item{alpha}{significance level}
  \item{n}{sample size}
  \item{del0}{upper limit set to \eqn{|p_{10}-p_{01}|} under the alternative hypothesis of
              equivalence, with \eqn{p_{10}} and \eqn{p_{01}} denoting the probabilities of
              discordant pairs of both kinds}
  \item{sw}{width of the search grid for determining the maximum of the rejection probability
           on the common boundary of the hypotheses}
  \item{tol}{upper bound to the absolute difference between size and target level below
             which the search for a corrected nominal level terminates}
  \item{maxh}{maximum number of interval halving steps to be carried out in finding the
              maximally raised nominal level}
}
\value{
  \item{alpha}{significance level}
  \item{n}{sample size}
  \item{del0}{upper limit set to \eqn{|p_{10}-p_{01}|} under the alternative hypothesis of
              equivalence, with \eqn{p_{10}} and \eqn{p_{01}} denoting the probabilities of
              discordant pairs of both kinds}
  \item{sw}{width of the search grid for determining the maximum of the rejection probability
           on the common boundary of the hypotheses}
   \item{ALPH_0}{value of the corrected nominal level obtained after nh steps}
   \item{NH}{number of interval-halving steps actually performed}
   \item{ERROR}{error indicator messaging "!!!!!" if the sufficient condition
                for the correctness of the result output by the program was found violated}
 }
\references{
Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 5.2.2.
}
\author{
Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>
}
\examples{
mcnemasc(0.05,50,0.20,0.05,0.0005,5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{paired samples}
\keyword{binary observations}
\keyword{McNemar setting}
\keyword{equivalence}
\keyword{asymptotic test}
\keyword{corrected nominal level}

