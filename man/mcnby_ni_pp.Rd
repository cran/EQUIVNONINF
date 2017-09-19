\name{mcnby_ni_pp}
\alias{mcnby_ni_pp}
\title{Computation of the posterior probability of the alternative hypothesis
of noninferiority in the McNemar setting, given a specific point in the sample space
}
\description{Evaluation of the integral on the right-hand side of Equation (5.24) 
on p. 88 of Wellek S (2010) Testing statistical hypotheses of equivalence and 
noninferiority. Second edition.
}
\usage{
mcnby_ni_pp(N,DEL0,N10,N01)
}
\arguments{
  \item{N}{sample size}
  \item{DEL0}{noninferiority margin to the difference of the parameters of the marginal
              binomial distributions under comparison}
  \item{N10}{count of pairs with \eqn{(X,Y) = (1,0)}}
  \item{N01}{count of pairs with \eqn{(X,Y) = (0,1)}}
}
\details{The program uses 96-point Gauss-Legendre quadrature on each of
10 subintervals into which the range of integration is partitioned.}
\value{
  \item{N}{sample size}
  \item{DEL0}{noninferiority margin to the difference of the parameters of the marginal
              binomial distributions under comparison}
  \item{N10}{count of pairs with \eqn{(X,Y) = (1,0)}}
  \item{N01}{count of pairs with \eqn{(X,Y) = (0,1)}}
\item{PPOST}{posterior probability of the alternative hypothesis \eqn{K_1: \delta > -\delta_0}
             with respect to the noninformative prior determined according to Jeffrey's rule}

}
\references{Wellek S: Testing statistical hypotheses of equivalence and 
        noninferiority. Second edition. Boca Raton: 
        Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 5.2.3.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\note{The program uses Equation (5.24) of Wellek S (2010) corrected for
a typo in the middle line which must read 
\deqn{\int_{\delta_0}^{(1+\delta_0)/2}\Big[ B\big(n_{01}+1/2,n-n_{01}+1\big)\,\,
p_{01}^{n_{01}-1/2}(1-p_{01})^{n-n_{01}}}.
}

\examples{
mcnby_ni_pp(72,0.05,4,5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{McNemar setting} 
\keyword{noninferiority}
\keyword{difference between response probabilities}
\keyword{posterior probability of the alternative hypothesis}
\keyword{Jeffrey's rule}
\keyword{Dirichlet prior}
\keyword{Gauss-Legendre quadrature}

