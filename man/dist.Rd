\name{dist}
\alias{dist}
\alias{dist.default}
\title{Distances in variouse approaches}
\description{
  Calculates a distance matrix from a dataset. 
}
\usage{
dist(x,\dots)
\method{dist}{default}(x,\dots)
          }
\arguments{
  \item{x}{a  dataset}
  \item{\dots}{further arguments to \code{\link[stats]{dist}}}
}
\value{
  a distance matrix
}
\details{
  The distance is computed based on \code{\link{cdt}}
  }
\seealso{
  \code{\link{aplus}}
}
\examples{
data(iris)
dist(iris[,1:4])
data(SimulatedAmounts)
dist(acomp(sa.lognormals),method="manhattan")
dist(rcomp(sa.lognormals))
dist(aplus(sa.lognormals))
dist(rplus(sa.lognormals))
}
\keyword{multivariate}
