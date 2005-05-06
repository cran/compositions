\name{totals}
\alias{totals}
\alias{totals.acomp}
\alias{totals.aplus}
\alias{totals.rcomp}
\alias{totals.rplus}
\title{Total sum of amounts }
\description{
  Calculates the total amount by summing the individual parts. 
}
\usage{
totals(x,\dots)
\method{totals}{acomp}(x,\dots)
\method{totals}{rcomp}(x,\dots)
\method{totals}{aplus}(x,\dots)
\method{totals}{rplus}(x,\dots)
          }
\arguments{
  \item{x}{an amount/amount dataset}
  \item{\dots}{not used, only here for generics}
}
\value{
  a numeric vector containing the total amounts
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{aplus}}
}
\examples{
data(SimulatedAmounts)
totals(acomp(sa.lognormals))
totals(rcomp(sa.lognormals,total=100))
totals(aplus(sa.lognormals))
totals(rplus(sa.lognormals))
aplus(acomp(sa.lognormals),total=totals(aplus(sa.lognormals)))
}
\keyword{multivariate}
