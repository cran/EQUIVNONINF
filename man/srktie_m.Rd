\name{srktie_m}
\alias{srktie_m}
\title{Analogue of srktie_d for settings where the distribution of intraindividual
differences is concentrated on a finite lattice}
\description{Analogue of the function srktie_d tailored for settings where the distribution
of the within-subject differences is concentrated on a finite lattice. For details see 
Wellek S (2010) Testing statistical hypotheses of equivalence and noninferiority. Second edition, pp.112-3.
}
\usage{
srktie_m(n,alpha,eps1,eps2,w,d) 
}
\arguments{
   \item{n}{sample size}
   \item{alpha}{significance level}
   \item{eps1}{absolute value of the left-hand limit of the hypothetical equivalence range for
                \eqn{q_+/(1-q_0) - 1/2}}
   \item{eps2}{right-hand limit of the hypothetical equivalence range for \eqn{q_+/(1-q_0) - 1/2}}
   \item{w}{span of the lattice in which the intraindividual differences take their values}
   \item{d}{row vector with the intraindividual differences for all \eqn{n} pairs as components}
 }
\details{Notation: \eqn{q_+} and \eqn{q_0} stands for the functional defined by 
         \eqn{q_+ = P[D_i+D_j>0]} and  \eqn{q_0 = P[D_i+D_j=0]}, respectively, 
         with \eqn{D_i} and \eqn{D_j} as the intraindividual differences observed in two individuals
         independently selected from the underlying bivariate population.
}
\value{
   \item{n}{sample size}
   \item{alpha}{significance level}
   \item{eps1}{absolute value of the left-hand limit of the hypothetical equivalence range for
                \eqn{q_+/(1-q_0) - 1/2}}
   \item{eps2}{right-hand limit of the hypothetical equivalence range for \eqn{q_+/(1-q_0) - 1/2}}
   \item{w}{span of the lattice in which the intraindividual differences take their values}
  \item{U_pl}{observed value of the \eqn{U}-statistics estimator of \eqn{q_+}}
  \item{U_0}{observed value of the \eqn{U}-statistics estimator of \eqn{q_0}}
  \item{UAS_PL}{observed value of \eqn{U_+/(1-U_0)}}
  \item{TAUHAS}{square root of the estimated asymtotic variance of \eqn{\sqrt{n}U_+/(1-U_0)}}
  \item{CRIT}{upper critical bound to \eqn{\sqrt{n}|U_+/(1-U_0) - 1/2 - 
               (\varepsilon_2-\varepsilon_1)/2|/\hat{\tau}}}
  \item{REJ}{indicator of a positive [=1] vs negative [=0] rejection decision to be taken with
             the data under analysis} 
  }
\references{
Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, pp. 112-114.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
d <- c(0.8,0.2,0.0,-0.1,-0.3,0.3,-0.1,0.4,0.6,0.2,0.0,-0.2,-0.3,0.0,0.1,0.3,-0.3,
       0.1,-0.2,-0.5,0.2,-0.1,0.2,-0.1)
srktie_m(24,0.05,0.2602,0.2602,0.1,d)
}
\keyword{Walsh averages}
\keyword{nonparametric paired-sample equivalence test}
\keyword{tied observations}
\keyword{finite lattice}
\keyword{U-statistics estimators}
\keyword{asymptotic normality}


