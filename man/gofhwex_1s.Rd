\name{gofhwex_1s}
\alias{gofhwex_1s}
\title{Critical constants of the exact UMPU test for absence of a substantial deficit of
heterozygotes as compared with a HWE-compliant SNP genotype distribution 
[noninferiority version of the test implemented by means of gofhwex]}
\description{The function computes the critical constants defining the UMPU test for 
one-sided equivalence of the population distribution of a SNP, to a distribution 
being in Hardy-Weinberg equilibrium (HWE).\cr
A substantial deficit of heterozygotes is defined to occur when the true value of the
parametric function \eqn{\omega = \frac{\pi_2/2}{\sqrt{\pi_1\pi_3}}} [called relative excess
heterozygosity (REH)] falls below unity by more than some given margin \eqn{\delta_0}.\cr
Like its two-sided counterpart [see the description of the R function gofhwex],
the test is conditional on the total count \eqn{S} of alleles of the kind of interest. 
}
\usage{
gofhwex_1s(alpha,n,s,del0) 
}
\arguments{
   \item{alpha}{significance level}
   \item{n}{number of genotyped individuals}
   \item{s}{observed count of alleles of the kind of interest}
   \item{del0}{noninferiority margin for \eqn{\omega}, which has to satisfy
         \eqn{\omega > 1-\delta_0} under the alternative hypothesis to be established}
 }
\value{
   \item{alpha}{significance level}
   \item{n}{number of genotyped individuals}
   \item{s}{observed count of alleles of the kind of interest}
   \item{del0}{noninferiority margin for \eqn{\omega}, which has to satisfy
         \eqn{\omega > 1-\delta_0} under the alternative hypothesis to be established}
  \item{C}{left-hand limit of the critical interval for the observed number \eqn{X_2} of heterozygotes} 
  \item{GAM}{probability of rejecting the null hypothesis when it turns out that \eqn{X_2=C}}
}
\references{
Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. 
Boca Raton: Chapman & Hall/CRC Press, 2010, pp. 300-302.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}

\examples{
gofhwex_1s(0.05,133,65,1-1/1.96)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{single nucleotide polymorphism (SNP)}
\keyword{Hardy-Weinberg equilibrium (HWE)}
\keyword{relative excess heterozygosity (REH)}
\keyword{noninferiority}
\keyword{uniformly most powerful unbiased test}

