\name{EFA.MRFA-package}
\alias{EFA.MRFA-package}
\title{
Dimensionality Assesment using Minimum Rank Factor Analysis (MRFA)
}
\description{
Package for performing Parallel Analysis using Minimum Rank Factor Analysis (MRFA) . It also include a function to perform the MRFA  only and another function to compute the Greater Lower Bound step for estimating the variables communalities.
}
\details{
For more information about the methods used in each function, please go to each main page.
}
\value{
\item{\code{\link{parallelMRFA}}}{Performs Parallel Analysis using Minimum Rank Factor Analysis (MRFA).}
\item{\code{\link{hullEFA}}}{Performs Hull analysis for assessing the number of factors to retain.}
\item{\code{\link{mrfa}}}{Performs Minimum Rank Factor Analysis (MRFA) procedure.}
\item{\code{\link{GreaterLowerBound}}}{Estimates the communalities of the variables from a factor model.}
}
\references{
Devlin, S. J., Gnanadesikan, R., & Kettenring, J. R. (1981). Robust estimation of dispersion matrices and principal components. Journal of the American Statistical Association, 76, 354-362. \doi{10.1080/01621459.1981.10477654}

Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011). The Hull Method for Selecting the Number of Common Factors. Multivariate Behavioral Research, 46(2), 340-364. \doi{10.1080/00273171.2011.564527}

ten Berge, J. M. F., & Kiers, H. A. L. (1991). A numerical approach to the approximate and the exact minimum rank of a covariance matrix. Psychometrika, 56(2), 309-315. \doi{10.1007/BF02294464}

Ten Berge, J.M.F., Snijders, T.A.B. & Zegers, F.E. (1981). Computational aspects of the greatest lower bound to reliability and constrained minimum trace factor analysis. Psychometrika, 46, 201-213.

Timmerman, M. E., & Lorenzo-Seva, U. (2011). Dimensionality assessment of ordered polytomous items with parallel analysis. Psychological Methods, 16(2), 209-220. \doi{10.1037/a0023353}
}
\author{
David Navarro-Gonzalez

Urbano Lorenzo-Seva
}
\examples{
## Example 1:

## perform a Parallel Analysis using an example Database with only 5 random data sets and
## using the 90th percentile of distribution of the random data
parallelMRFA(IDAQ, Ndatsets=5, percent=90)

## For speeding purposes, the number of datasets have been largely reduced. For a proper
## use of parallelMRFA, we recommend to use the default Ndatsets value (Ndatsets=500)

#Example 2:

## Perform the Hull method defining the maximum number of dimensions to be tested by the
## Parallel Analysis + 1 rule, with Maximum Likelihood factor extraction method and CAF
## as Hull index.
hullEFA(IDAQ, extr = "ML")
}
\keyword{ package }
