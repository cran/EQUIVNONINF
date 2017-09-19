\name{EQUIVNONINF-package}
\alias{EQUIVNONINF-package}
\alias{EQUIVNONINF}
\docType{package}
\title{Testing for equivalence and noninferiority
}
\description{
The package makes available in R the complete set of programs accompanying S. Wellek's (2010) monograph
             "Testing Statistical Hypotheses of Equivalence and Noninferiority. Second Edition"
             (Chapman&Hall/CRC). 
}
\author{Stefan Wellek <stefan.wellek@zi-mannheim.de> \cr
        Peter Ziegler <peter.ziegler@zi-mannheim.de> \cr \cr
Maintainer: Stefan Wellek <stefan.wellek@zi-mannheim.de>
}
\references{Wellek S: Testing statistical hypotheses of equivalence and noninferiority. Second edition. Boca Raton: Chapman & Hall/CRC Press, 2015.

}

\keyword{ package }
\examples{
bi2ste1(397,397,0.0,0.025,0.511,0.384)
bi2ste2(0.0,0.025,0.95,0.8,0.80,1.0)
}
\note{In order to keep execution time of all examples below the limit set by the
CRAN administration, in a number of cases the function calls shown in the documentation contain 
specifications which are insufficient for real applications. This holds in particular true for the 
width sw of search grids, which should be chosen to be .001 or smaller. 
Similarly, the maximum number of interval halving steps to be carried out in finding maximally 
admissible significance levels should be set to values >= 10.}

