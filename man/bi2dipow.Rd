\name{bi2dipow}
\alias{bi2dipow}
\title{Exact rejection probability of the asymptotic test for equivalence of two
unrelated binomial proportions with respect to the difference of their expectations
at any nominal level under an arbitrary parameter configuration
}
\description{The program computes exact values of the rejection probability of the asymptotic
test for equivalence in the sense of  \eqn{-\delta_1 < p_1-p_2 < \delta_2}, at any nominal
level \eqn{\alpha_0}. [The largest \eqn{\alpha_0} for which the test is valid in terms of the
significance level, can be computed by means of the program bi2diffac.]
}
\usage{
bi2dipow(alpha0,m,n,del1,del2,p1,p2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha0}{nominal significance level}
  \item{m}{size of Sample 1}
  \item{n}{size of Sample 2}
  \item{del1}{absolute value of the lower limit of the hypothetical equivalence range for \eqn{p_1-p_2}}
  \item{del2}{upper limit of the hypothetical equivalence range for \eqn{p_1-p_2}}
  \item{p1}{true value of the success probability in Population 1}
  \item{p2}{true value of the success probability in Population 2}

}

\value{
  \item{alpha0}{nominal significance level}
  \item{m}{size of Sample 1}
  \item{n}{size of Sample 2}
  \item{del1}{absolute value of the lower limit of the hypothetical equivalence range for \eqn{p_1-p_2}}
  \item{del2}{upper limit of the hypothetical equivalence range for \eqn{p_1-p_2}}
  \item{p1}{true value of the success probability in Population 1}
  \item{p2}{true value of the success probability in Population 2}
   \item{POWEX0}{exact rejection probability under \eqn{(p_1,p_2)} of the test at nominal level \eqn{\alpha_0}
                 for equivalence of two binomial distributions with respect to the difference of
                 the success probabilities}
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
bi2dipow(0.0228,50,50,0.20,0.20,0.50,0.50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{binomial two-sample problem}
\keyword{equivalence}
\keyword{difference of success probabilities}
\keyword{asymptotic test}
\keyword{exact power}

