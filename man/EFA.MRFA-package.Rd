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
\item{\code{\link{testme}}}{An auto-executable script for testing the functions included in DA.MRFA.}
}
\references{
Devlin, S. J., Gnanadesikan, R., & Kettenring, J. R. (1981). Robust estimation of dispersion matrices and principal components. Journal of the American Statistical Association, 76, 354-362. \url{http://doi.org/10.1080/01621459.1981.10477654}

Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011). The Hull Method for Selecting the Number of Common Factors. Multivariate Behavioral Research, 46(2), 340-364. \url{https://doi.org/10.1080/00273171.2011.564527}

ten Berge, J. M. F., & Kiers, H. A. L. (1991). A numerical approach to the approximate and the exact minimum rank of a covariance matrix. Psychometrika, 56(2), 309–315. \url{http://doi.org/10.1007/BF02294464}

Ten Berge, J.M.F., Snijders, T.A.B. & Zegers, F.E. (1981). Computational aspects of the greatest lower bound to reliability and constrained minimum trace factor analysis. Psychometrika, 46, 201-213.

Timmerman, M. E., & Lorenzo-Seva, U. (2011). Dimensionality assessment of ordered polytomous items with parallel analysis. Psychological Methods, 16(2), 209-220. \url{http://doi.org/10.1037/a0023353}
}
\author{
David Navarro-Gonzalez

Urbano Lorenzo-Seva
}
\examples{
## Each man page contains examples of each function. For a fast global example use
testme(example = TRUE)

## For speeding purposes, the number of datasets have been largely reduced. For a proper
## use of testme, use it without the "example" argument or using example=FALSE.
}
\keyword{ package }
