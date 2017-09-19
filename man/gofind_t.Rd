\name{gofind_t}
\alias{gofind_t}
\title{Establishing approximate independence in a two-way contingency table: Test statistic
and critical bound}
\description{The function computes all quantities required for carrying out the asymptotic test
for approximate independence of two categorial variables derived in \eqn{\S}{Par.} 9.2 of
Wellek S (2010) Testing statistical hypotheses of equivalence and noninferiority. Second edition.
}
\usage{
gofind_t(alpha,r,s,eps,xv) 
}
\arguments{
   \item{alpha}{significance level}
   \item{r}{number of rows    of the contingency table under analysis}
   \item{s}{number of columns of the contingency table under analysis}
   \item{eps}{margin to the Euclidean distance between the vector \eqn{\mathbf{\pi}} of true cell 
              probabilities and the associated vector of products of marginal totals}
  \item{xv}{row vector of length \eqn{r\cdot s}{r * s} whose \eqn{(i-1)s + j}-th component is the entry
            in cell \eqn{(i,j)} of the \eqn{r\times s}{r * s} contingency table under analysis \eqn{i=1,\ldots,r},
            \eqn{j=1,\ldots,s}.}
}
\value{
  \item{n}{size of the sample to which the input table relates}
   \item{alpha}{significance level}
   \item{r}{number of rows    of the contingency table under analysis}
   \item{s}{number of columns of the contingency table under analysis}
   \item{eps}{margin to the Euclidean distance between the vector \eqn{\mathbf{\pi}} of true cell 
              probabilities and the associated vector of products of marginal totals}
  \item{X(r,s)}{observed cell counts}
  \item{DSQ_OBS}{observed value of the squared Euclidean distance}
  \item{VN}{square root of the estimated asymtotic variance of \eqn{\sqrt{n}DSQ\_OBS}}
  \item{CRIT}{upper critical bound to \eqn{\sqrt{n}DSQ\_OBS}}
  \item{REJ}{indicator of a positive [=1] vs negative [=0] rejection decision to be taken with
             the data under analysis}
}
\references{
Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 9.2.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
xv <- c(8, 13, 15, 6,  19, 21, 31, 7) 
gofind_t(0.05,2,4,0.15,xv)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{2x2 contingency table}
\keyword{establishing approximate independence of two categorical variables}
\keyword{Euclidean distance statistic}
\keyword{asymptotic normality}


