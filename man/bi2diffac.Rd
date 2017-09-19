\name{bi2diffac}
\alias{bi2diffac}
\title{Determination of a corrected nominal significance level for the asymptotic test for equivalence
of two unrelated binomial proportions with respect to the difference \eqn{\delta} 
of their population counterparts
}
\description{The program computes the largest nominal significance level 
which can be substituted for the target level \eqn{\alpha} without making the exact 
size of the asymptotic testing procedure larger than \eqn{\alpha}.  } 
\usage{
bi2diffac(alpha,m,n,del1,del2,sw,tolrd,tol,maxh)
}
\arguments{
  \item{alpha}{significance level}
  \item{m}{size of Sample 1}
  \item{n}{size of Sample 2}
  \item{del1}{absolute value of the lower limit of the hypothetical equivalence range for \eqn{p_1-p_2}}
  \item{del2}{upper limit of the hypothetical equivalence range for \eqn{p_1-p_2}}
  \item{sw}{width of the search grid for determining the maximum of the rejection probability
           on the common boundary of the hypotheses}
  \item{tolrd}{horizontal distance of the left- and right-most boundary point to be included 
               in the search grid}
  \item{tol}{upper bound to the absolute difference between size and target level below
             which the search for a corrected nominal level terminates}
  \item{maxh}{maximum number of interval halving steps to be carried out in finding the
              maximally raised nominal level}
}
\value{
  \item{alpha}{significance level}
  \item{m}{size of Sample 1}
  \item{n}{size of Sample 2}
  \item{del1}{absolute value of the lower limit of the hypothetical equivalence range for \eqn{p_1-p_2}}
  \item{del2}{upper limit of the hypothetical equivalence range for \eqn{p_1-p_2}}
  \item{sw}{width of the search grid for determining the maximum of the rejection probability
           on the common boundary of the hypotheses}
  \item{tolrd}{horizontal distance of the left- and right-most boundary point to be included 
               in the search grid}
  \item{tol}{upper bound to the absolute difference between size and target level below
             which the search for a corrected nominal level terminates}
  \item{maxh}{maximum number of interval halving steps to be carried out in finding the
              maximally raised nominal level}
   \item{NH}{number of interval-halving steps actually performed}
   \item{ALPH_0}{value of the raised nominal level obtained after NH steps}
   \item{SIZE0}{size of the critical region corresponding to \eqn{\alpha_0}{alpha_0}}
   \item{ERROR}{error indicator answering the question of whether or not the sufficient condition
                for the correctness of the result output by the program, was satisfied} }
\references{
Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 6.6.6.
}
\author{
Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>
}
\examples{
bi2diffac(0.05,20,20,0.40,0.40,0.1,1e-6,1e-4,3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{binomial two-sample problem}
\keyword{equivalence}
\keyword{difference of success probabilities}
\keyword{asymptotic test}
\keyword{corrected nominal level}

