\name{mcnempow}
\alias{mcnempow}
\title{Exact rejection probability of the asymptotic test for equivalence of two
paired binomial proportions with respect to the difference of their expectations 
(McNemar setting)}
\description{The program computes exact values of the rejection probability of the asymptotic
test for equivalence in the sense of  \eqn{-\delta_0 < p_{10}-p_{01} < \delta_0}, at any nominal
level \eqn{\alpha}. [The largest \eqn{\alpha} for which the test is valid in terms of the
significance level, can be computed by means of the program mcnemasc.] }
\usage{
 mcnempow(alpha,n,del0,p10,p01)
}
\arguments{
  \item{alpha}{nominal significance level}
  \item{n}{sample size}
    \item{del0}{upper limit set to \eqn{|\delta|} under the alternative hypothesis of
              equivalence}
  \item{p10}{true value of \eqn{P[X=1,Y=0]}}
  \item{p01}{true value of \eqn{P[X=0,Y=1]}}
}
\value{
  \item{alpha}{nominal significance level}
  \item{n}{sample size}
    \item{del0}{upper limit set to \eqn{|\delta|} under the alternative hypothesis of
              equivalence}
  \item{p10}{true value of \eqn{P[X=1,Y=0]}}
  \item{p01}{true value of \eqn{P[X=0,Y=1]}}
   \item{POW}{exact rejection probability of the asymptotic McNemar test for equivalence 
              at nominal level \eqn{\alpha}}
    \item{ERROR}{error indicator messaging "!!!!!" if the sufficient condition
                for the correctness of the result output by the program was found violated}
   }
\references{
Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, p.84.
}
\author{
Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>
}
\examples{
mcnempow(0.024902,50,0.20,0.30,0.30)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{McNemar setting}
\keyword{asymptotic test for equivalence}
\keyword{exact power}

