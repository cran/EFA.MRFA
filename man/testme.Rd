\name{testme}
\alias{testme}
\title{
Function for test the full capabilites of \code{EFA.MRFA} package
}
\description{
testme is a function which includes all the functions available in the EFA.MRFA package, and his purpose is to serve as a demonstration of the potential usage of the package functions.
}
\usage{
testme(example = FALSE)
}
\arguments{
  \item{example}{
  Determine if the usage is for exampling purpose or not (FALSE by default).
}
}
\details{
The function implements all the functions on the \code{EFA.MRFA} package. For ilustrating purposes, the loading matrix obtained using MRFA (\code{\link{mrfa}} function) is rotated using the function \code{\link[PCovR]{promin}} (Lorenzo-Seva, 1999), from \code{PCovR} package (Vervloet, Kiers, Van den Noortgate & Ceulemans, 2015).

Both procedures for assessing the number of factors to be retained are provided. As can be seen, the procedures differ on the number of factors suggested, with Parallel Analysis suggesting 3 factors (corresponding to the 3 sub-dimensions of aggression), while Hull only suggests one factor to be retained, corresponding to the general agression factor.
}
\references{
Lorenzo-Seva, U. (1999). Promin: A Method for Oblique Factor Rotation. Multivariate Behavioral Research, 34(3), 347–365 \url{http://doi.org/10.1207/S15327906MBR3403_3}

Vervloet, M., Kiers, H. A. L., Van den Noortgate, W., & Ceulemans, E. (2015). PCovR: An R Package for Principal Covariates Regression. JSS Journal of Statistical Software, 65(8). \url{http://doi.org/10.18637/jss.v065.i08}
}
\author{
David Navarro-Gonzalez

Urbano Lorenzo-Seva
}
\examples{
## For a fast usage of testme use:
testme(example=TRUE)

## For speeding purposes, the number of datasets have been largely reduced. For a proper
## use of testme, use testme() or using testme(example=FALSE).

}
