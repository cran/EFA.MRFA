\name{IDAQ}
\alias{IDAQ}
\docType{data}
\title{
IDAQ database
}
\description{
A database to be used as example in the functions included on \code{DA.MRFA} package. It contains the answers of 100 participants to IDAQ questionnaire (Ruiz-Pamies, Lorenzo-Seva, Morales-Vives, Cosi, Vigil-Colet, 2014), which was developed for assessing Physical, Verbal and Indirect aggression. The original questionnaire contains 27 Likert-items, ranging from 1 to 5.
}
\usage{data("IDAQ")}
\format{
  A data frame with 100 observations and 23 variables measuring 3 different types of aggression (Physical, Verbal and Indirect).
}
\details{
The original sample contains 27 items, because includes 4 Social Desirability markers, but for the purpose of the \code{DA.MRFA} functions, they had been removed. Also, the original sample contains 750 participants, and the following database only contains 100 for speeding purposes.
}
\source{
More information about the questionnaire can be found at:

\url{http://psico.fcep.urv.cat/tests/idaq/en/descripcion.html}
}
\references{
Ruiz-Pamies, M., Lorenzo-Seva, U., Morales-Vives, F., Cosi, S., & Vigil-Colet, A. (2014). I-DAQ: a new test to assess direct and indirect aggression free of response bias. The Spanish Journal of Psychology, 17, E41. \doi{10.1017/sjp.2014.43}
}
\examples{
data(IDAQ)
}
\keyword{datasets}
