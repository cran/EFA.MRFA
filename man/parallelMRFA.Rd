\name{parallelMRFA}
\alias{parallelMRFA}

\title{
Parallel Analysis using Minimum Rank Factor Analysis (MRFA)
}
\description{
Performs Parallel Analysis using Minimum Rank Factor Analysis (MRFA).
}
\usage{
parallelMRFA(X, Ndatsets = 500, percent = 95, corr= "Pearson", display = TRUE,
    graph = TRUE)
}
\arguments{
  \item{X}{
Raw sample scores.
}
  \item{Ndatsets}{
Number of random datasets used to compute the random distribution of eigenvalues.
}
  \item{percent}{
Desired percentile of distribution of random eigenvalues (for example 95 for the 95th percentile) to be used as threshold.
}
  \item{corr}{
Determine if Pearson or Polychoric matrix will be used "Pearson": Computes Pearson correlation matrix "Polychoric": Computes Polychoric/Tetrachoric correlation matrix (heavy time consuming).
}
  \item{display}{
Determines if the output will be displayed in the console, TRUE by default. If it is TRUE, the output is returned silently and if it is FALSE, the output is returned in the console.
}
  \item{graph}{
Request a plot representing the percentage of explained variance by the real data, by the mean of the random data and using the percentile of distribution of random eigenvalues, defined in the percent argument.
}
}
\details{
\code{parallelMRFA} is based on the procedure proposed by Timmerman and Lorenzo-Seva (2011) which is designed for assessing the dimensionality of a variable set. The principal advantage of using MRFA (Ten Berge & Kiers, 1991) instead the usual PCA extraction process is that the eigenvalues obtained from MRFA can be used to estimate the explained common variance per factor.

The eigenvalue sampling distribution is obtaining using a nonparametric approach: a permutation of the raw data (Buja & Eyuboglu, 1992). This approach is recommended for PA especially in cases where the observed data ditribution clearly deviates from normality.

If the matrix to analyze is not positive-defined, a smoothering procedure will be applied (Devlin, Gnanadesikan & Kettenring, 1981).
}
\value{
\item{Real_Data}{A vector containing the percentage of explained variance by the real data for each factor}
\item{Mean_random}{A vector containing the percentage of explained variance by the mean of random data for each factor}
\item{Percentile_random}{A vector containing the percentage of explained variance by the percentile of distribution of random data for each factor}
\item{Number_factors_mean}{The number of factors to be retained suggested comparing the real data with the mean of the random data}
\item{Number_factors_percentiles}{The number of factors to be retained suggested comparing the real data with the percentile of distribution of the random data}
}
\references{
Buja, A., & Eyuboglu, N. (1992). Remarks on Parallel Analysis. Multivariate Behavioral Research, 27(4), 509-540. \doi{10.1207/s15327906mbr2704_2}

Devlin, S. J., Gnanadesikan, R., & Kettenring, J. R. (1981). Robust estimation of dispersion matrices and principal components. Journal of the American Statistical Association, 76, 354-362. \doi{10.1080/01621459.1981.10477654}

ten Berge, J. M. F., & Kiers, H. A. L. (1991). A numerical approach to the approximate and the exact minimum rank of a covariance matrix. Psychometrika, 56(2), 309–315. \doi{10.1007/BF02294464}

Timmerman, M. E., & Lorenzo-Seva, U. (2011). Dimensionality assessment of ordered polytomous items with parallel analysis. Psychological Methods, 16(2), 209-220. \doi{10.1037/a0023353}
}
\author{
David Navarro-Gonzalez

Urbano Lorenzo-Seva
}
\examples{
## perform a Parallel Analysis using an example Database with only 10 random data sets and
## using the 90th percentile of distribution of the random data
parallelMRFA(IDAQ, Ndatsets=10, percent=90)

## For speeding purposes, the number of datasets have been largely reduced. For a proper
## use of parallelMRFA, we recommend to use the default Ndatsets value (Ndatsets=500)

}
