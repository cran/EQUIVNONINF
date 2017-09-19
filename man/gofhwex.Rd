\name{gofhwex}
\alias{gofhwex}
\title{Critical constants of the exact UMPU test for approximate compatibility of a
SNP genotype distribution with the Hardy-Weinberg model}
\description{The function computes the critical constants defining the uniformly most 
powerful unbiased test for equivalence of the population distribution of the three genotypes 
distinguishable in terms of a single nucleotide polymorphism (SNP), to a distribution 
being in Hardy-Weinberg equilibrium (HWE).\cr
The test is conditional on the total count \eqn{S} of alleles of the kind of interest, and
the parameter \eqn{\theta}, in terms of which equivalence shall be established, is defined
by \eqn{\theta = \frac{\pi_2^2}{\pi_1(1-\pi_1-\pi_2)}}, with \eqn{\pi_1} and \eqn{\pi_2} denoting
the population frequence of homozygotes of the 1st kind and heterozygotes, respectively.
}
\usage{
gofhwex(alpha,n,s,del1,del2) 
}
\arguments{
   \item{alpha}{significance level}
   \item{n}{number of genotyped individuals}
   \item{s}{observed count of alleles of the kind of interest}
   \item{del1}{absolute value of the lower equivalence limit to \eqn{\theta/4 - 1}}
   \item{del2}{upper equivalence limit to \eqn{\theta/4 - 1}}
}
\value{
   \item{alpha}{significance level}
   \item{n}{number of genotyped individuals}
   \item{s}{observed count of alleles of the kind of interest}
   \item{del1}{absolute value of the lower equivalence limit to \eqn{\theta/4 - 1}}
   \item{del2}{upper equivalence limit to \eqn{\theta/4 - 1}}
\item{C1}{left-hand limit of the critical interval for the observed number \eqn{X_2} of heterozygotes} 
\item{C2}{right-hand limit of the critical interval for the observed number \eqn{X_2}}
\item{GAM1}{probability of rejecting the null hypothesis when it turns out that \eqn{X_2=C_1}}
\item{GAM2}{probability of rejecting the null hypothesis for \eqn{X_2=C_2}}
}
\references{Wellek S: Tests for establishing compatibility of an observed genotype
distribution with Hardy-Weinberg equilibrium in the case of a biallelic
locus. Biometrics 60 (2004), 694-703.

Goddard KAB, Ziegler A, Wellek S: Adapting the logical basis of tests
for Hardy-Weinberg equilibrium to the real needs of association studies in human and
medical genetics. Genetic Epidemiology 33 (2009), 569-580.

Wellek S: Testing statistical hypotheses of equivalence and noninferiority. 
Second edition. Boca Raton: Chapman & Hall/CRC Press, 2010, \eqn{\S}{Par.} 9.4.2.}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
Peter Ziegler <peter.ziegler@zi-mannheim.de>}
\examples{
gofhwex(0.05,475,429,1-1/1.96,0.96)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{single nucleotide polymorphism (SNP)}
\keyword{Hardy-Weinberg equilibrium (HWE)}
\keyword{model validation}
\keyword{uniformly most powerful unbiased test}

