\name{gsiinternal}
\alias{gsi}
\title{Environment containing the the old gsi functions}
\description{
  If you need some of the old depricated gsi.functions you can attach
  this environment, which contains them all.
}
\usage{
gsi
}
\arguments{
  None
}
\value{
unkown
}
\details{
  Will not be given.
  
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\seealso{
  \code{\link{gsi}}
}
\examples{
attach(gsi)
objects("gsi")
detach(gsi)
"Hallo dont call these functions, they will be removed in the next version"
}
\keyword{internal}
