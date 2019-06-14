\name{GreaterLowerBound}
\alias{GreaterLowerBound}
\title{
Greater Lower Bound step (glb)
}
\description{
Estimates the communalities of the variables from a factor model where the number of factors is the number with positive eigenvalues.
}
\usage{
GreaterLowerBound(C, conv = 0.000001, T, pwarnings = FALSE)
}
\arguments{
  \item{C}{
Covariance/correlation matrix to be used in the analysis.
}
  \item{conv}{
Convergence criterion for glb step. The default convergence criterion will be conv=0.000001 . If the user determine a specific value, this will prevail.
}
  \item{T}{
Random matrix for start (can be omitted). If provided, it has to be the same size than the matrix provided in the C argument.
}
  \item{pwarnings}{
Determines if the possible warnings occurred during the computation will be printed in the console.
}
}
\details{
Code adapted from a MATLAB function by Jos Ten Berge based on Ten Berge, Snijders & Zegers (1981) and Ten Berge & Kiers (1991).
}
\value{
  \item{gam}{Optimal communalities for each variable}
}
\references{
Ten Berge, J.M.F., & Kiers, H.A.L. (1991). A numerical approach to the exact and the approximate minimum rank of a covariance matrix. Psychometrika, 56, 309-315.

Ten Berge, J.M.F., Snijders, T.A.B. & Zegers, F.E. (1981). Computational aspects of the greatest lower bound to reliability and constrained minimum trace factor analysis. Psychometrika, 46, 201-213.
}
\author{
David Navarro-Gonzalez

Urbano Lorenzo-Seva
}

\examples{
## perform glb using the correlation matrix of the IDAQ dataset, and using severe convergence
## criterion.
GreaterLowerBound(cor(IDAQ), conv=0.000001)
}
