\name{summary.aplus}
\alias{summary.aplus}
\alias{summary.rplus}
\alias{summary.rmult}
\title{Summaries of amounts }
\description{
  Summary of a vector of amounts, according to its underlying geometry.
}
\usage{
          \method{summary}{aplus}( object, \dots ,digits=max(3, getOption("digits")-3))
          \method{summary}{rplus}( object, \dots )
          \method{summary}{rmult}( object, \dots )
          }
\arguments{
  \item{object}{an \code{\link{aplus}}/\code{\link{rplus}} set of
    amounts}
  \item{digits}{the number of significant digits to be used. The
    argument can also be used with rplus/rmult.}
  \item{\dots}{not used, only here for generics}
}
\value{
  A matrix containing summary statistics (minimum, the three quantiles, the mean 
  and the maximum) of each component.
  The obtained value is the same as for the classical summary \code{\link[base]{summary}},
  although in the case of \code{\link{aplus}} objects, the statistics have been computed in a
  logarithmic geometry, and exponentiated afterwards (which just changes the mean, equivalent 
  to the geometric mean of the data set).
}
\details{
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{aplus}},\code{\link{rplus}},\code{\link{summary.acomp}},
  \code{\link{summary.rcomp}}
}
\references{
  

}
\examples{
data(SimulatedAmounts)
summary(aplus(sa.lognormals))
summary(rplus(sa.lognormals))
summary(rmult(sa.lognormals))

}
\keyword{multivariate}
