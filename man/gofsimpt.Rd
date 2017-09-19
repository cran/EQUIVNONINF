\name{gofsimpt}
\alias{gofsimpt}
\title{Establishing goodness of fit of an observed to a fully specified multinomial distribution:
test statistic and critical bound}
\description{The function computes all quantities required for carrying out the asymptotic test for goodness 
rather than lack of fit of an observed to a fully specified multinomial distribution
derived in \eqn{\S}{Par.} 9.1 of Wellek S (2010) Testing statistical hypotheses of equivalence and noninferiority. 
Second edition.
}
\usage{
gofsimpt(alpha,n,k,eps,x,pio) 
}
\arguments{
   \item{alpha}{significance level}
   \item{n}{sample size} 
   \item{k}{number of categories}
   \item{eps}{margin to the Euclidean distance between the vectors \eqn{\mathbf{\pi}} and
              \eqn{\mathbf{\pi}_0} of true and hypothesized cell probabilities}
   \item{x}{vector of length \eqn{k} with the observed cell counts as components} 
   \item{pio}{prespecified vector of cell probabilities}    
 }
\value{
   \item{alpha}{significance level}
   \item{n}{sample size} 
   \item{k}{number of categories}
   \item{eps}{margin to the Euclidean distance between the vectors \eqn{\mathbf{\pi}} and
              \eqn{\mathbf{\pi}_0} of true and hypothesized cell probabilities}
  \item{X(1,K)}{observed cell counts}
  \item{PI0(1,K)}{hypothecized cell probabilities}
  \item{DSQPIH_0}{observed value of the squared Euclidean distance}
  \item{VN_N}{square root of the estimated asymtotic variance of \eqn{\sqrt{n}DSQPIH\_0}}
  \item{CRIT}{upper critical bound to \eqn{\sqrt{n}DSQPIH\_0}}
  \item{REJ}{indicator of a positive [=1] vs negative [=0] rejection decision to be taken with
             the data under analysis}
 }
\references{
Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 9.1.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
x<- c(17,16,25,9,16,17)
pio <- rep(1,6)/6
gofsimpt(0.05,100,6,0.15,x,pio)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{one-dimensional contingency table}
\keyword{fully specified multinomial distribution}
\keyword{model verification}
\keyword{Euclidean distance statistic}
\keyword{asymptotic normality}


