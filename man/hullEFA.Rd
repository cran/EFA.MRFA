\name{hullEFA}
\alias{hullEFA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hull method for selecting the number of common factors
}
\description{
Performs the Hull method (Lorenzo-Seva, Timmerman, & Kiers, 2011), which aims to find a model with an optimal balance between model fit and number of parameters.
}
\usage{
hullEFA(X, maxQ, extr = "ULS", index_hull = "CAF", display = TRUE, graph = TRUE,
        details = TRUE)
}
\arguments{
  \item{X}{
Raw sample scores.
}
  \item{maxQ}{
Maximum of dimensions to be tested. By default it will be determined by the Parallel Analysis advised dimensions, but the user can define it manually.
}
  \item{extr}{
Extraction method, the two options available being: "ULS" (Unweigthed Least Squares, by default) and "ML" (Maximum Likelihood).
}
  \item{index_hull}{
The index that will be used for determining the number of dimensions. The available options are the following: "CAF", "CFI""RMSEA", being "CAF" by default.
}
  \item{display}{
Determines if the output will be displayed in the console, TRUE by default. If it is TRUE, the output is returned silently and if it is FALSE, the output is returned in the console.
}
  \item{graph}{
Request a plot representing the Hull curve.
}
  \item{details}{
If detailed table will be displayed, containing the factors outside the convex Hull.
  }
}
\details{
\code{hullEFA} is based on the procedure proposed by Lorenzo-Seva, Timmerman, &  Kiers (2011) which is designed for assessing the dimensionality of a variable set. The hull heuristic was originally proposed by Ceulemans & Kiers (2006) in the context of model selection in multiway data analysis.

The hull analysis is performed in four main steps:

1. The range of factors to be considered is determined.

2. The goodness-of-fit of a series of factor solutions is assessed.

3. The degrees of freedom of the series of factor solutions is computed.

4. The elbow is located in the higher boundary of the convex hull of the hull plot.


The number of factors extracted in the solution associated with the elbow is considered the optimal number of common factors.

In the Lorenzo-Seva, Timmerman, & Kiers (2011) simulation study, the Hull method outperformed the other selected methods in recovering the corrrect number of major factors.
}
\value{
  \item{Matrix}{Matrix containing the results of the Hull method using for the selected index.}
  \item{n_factors}{Number of advised dimensions by the selected index.}
}
\references{
Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011). The Hull Method for Selecting the Number of Common Factors. Multivariate Behavioral Research, 46(2), 340-364. \doi{10.1080/00273171.2011.564527}

Ceulemans, E., & Kiers, H. A. L. (2006). Selecting among three-mode principal component models of different types and complexities: A numerical convex hull based method. British Journal of Mathematical and Statistical Psychology, 59: 133–150. \doi{10.1348/000711005X64817}
}
\author{
David Navarro-Gonzalez

Urbano Lorenzo-Seva
}
\examples{
## Perform the Hull method defining the maximum number of dimensions to be tested by the
## Parallel Analysis + 1 rule, with Maximum Likelihood factor extraction method and CAF
## as Hull index.
hullEFA(IDAQ, extr = "ML")
}
