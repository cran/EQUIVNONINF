\name{bi2st}
\alias{bi2st}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Critical constants for the exact Fisher type UMPU test for equivalence of two binomial distributions with respect to
the odds ratio}
\description{The function computes the critical constants defining the uniformly most powerful unbiased test for
equivalence of two binomial distributions with parameters
\eqn{p_1} and \eqn{p_2} in terms of the odds ratio. 
Like the ordinary Fisher type test of the null hypothesis
\eqn{p_1 = p_2}, the test is conditional on the total number 
\eqn{S} of successes in the pooled sample.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bi2st(alpha,m,n,s,rho1,rho2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{alpha}{significance level}
   \item{m}{size of Sample 1}
   \item{n}{size of Sample 2}
   \item{s}{observed total count of successes}
   \item{rho1}{lower limit of the hypothetical equivalence
            range for the odds ratio 
            \eqn{\varrho = \frac{p_1(1-p_2)}{p_2(1-p_1)}}}
   \item{rho2}{upper limit of the hypothetical equivalence
            range for \eqn{\varrho}}
}
\value{
   \item{alpha}{significance level}
   \item{m}{size of Sample 1}
   \item{n}{size of Sample 2}
   \item{s}{observed total count of successes}
   \item{rho1}{lower limit of the hypothetical equivalence
            range for the odds ratio 
            \eqn{\varrho = \frac{p_1(1-p_2)}{p_2(1-p_1)}}}
   \item{rho2}{upper limit of the hypothetical equivalence
            range for \eqn{\varrho}}
\item{C1}{left-hand limit of the critical interval for
          the number \eqn{X} of successes observed in
          Sample 1}
\item{C2}{right-hand limit of the critical interval for
          \eqn{X}}
\item{GAM1}{probability of rejecting the null hypothesis 
            when it turns out that \eqn{X=C_1}}
\item{GAM2}{probability of rejecting the null hypothesis 
            for \eqn{X=C_2}}
}
\references{Wellek S: Testing statistical hypotheses of equivalence and noninferiority. 
Second edition. Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 6.6.4.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
bi2st(.05,225,119,171, 2/3, 3/2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{binomial two-sample problem}
\keyword{equivalence}
\keyword{odds ratio}
\keyword{exact Fisher-type test}
\keyword{UMPU test}
\keyword{conditional test}
\keyword{randomized test}
