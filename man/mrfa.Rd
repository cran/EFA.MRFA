\name{mrfa}
\alias{mrfa}
\title{
Minimum Rank Factor Analysis function
}
\description{
Performs Minimum Rank Factor Analysis (MRFA) procedure, proposed by Ten Berge & Kiers (1991).
}
\usage{
mrfa(SIGMA, dimensionality = 1, random = 10, conv1, conv2, display = TRUE,
    pwarnings = FALSE)
}
\arguments{
  \item{SIGMA}{
Covariance/correlation matrix to be used in the analysis.
}
  \item{dimensionality}{
Common factors used to find communality estimates. The value has to be between 0 and the number of items minus 1, being the default option: 1 dimension to be retained. If 0 is selected, a more strict convergence criterion will be used.
}
  \item{random}{
Number of random starts.
}
  \item{conv1}{
Convergence criterion for MRFA. The default convergence criterion will be conv1=0.0001 . If the user determine a specific value, this will prevail.
}
  \item{conv2}{
Convergence criterion for glb step. The default convergence criterion will be conv2=0.001 . If the user determine a specific value, this will prevail.
}
  \item{display}{
Determines if the output will be displayed in the console, TRUE by default. If it is TRUE, the output is returned silently and if it is FALSE, the output is returned in the console.
}
  \item{pwarnings}{
Determines if the possible warnings occurred during the computation will be printed in the console.
}
}

\value{
  \item{A}{Factor loading matrix}
  \item{Matrix}{Covariance/Correlation matrix with optimal communalities in the diagonal}
  \item{gam}{Optimal communalities for each variable}
}
\references{
ten Berge, J. M. F., & Kiers, H. A. L. (1991). A numerical approach to the approximate and the exact minimum rank of a covariance matrix. Psychometrika, 56(2), 309-315. \doi{10.1007/BF02294464}
}
\author{
David Navarro-Gonzalez

Urbano Lorenzo-Seva
}
\examples{
## perform MRFA using the correlation matrix of the IDAQ dataset, and using the default
## convergence criterion for MRFA and glb step.
mrfa(cor(IDAQ), dimensionality=3)
}
