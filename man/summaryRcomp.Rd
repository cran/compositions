\name{summary.rcomp}
\alias{summary.rcomp}
\title{Summary of compositions in real geometry}
\description{
  Compute a summary of a composition based on real geometry.
}
\usage{
          \method{summary}{rcomp}( object, \dots )
          }
\arguments{
  \item{object}{an \code{\link{rcomp}} dataset of compositions }
  \item{\dots}{further arguments to \code{summary}}
}
\value{
  A matrix containing summary statistics.
  The value is the same as for the classical summary
  \code{\link[base]{summary}} applied to  a closed dataset.

}
\details{
  The data is applied a \code{\link{clo}} operation before the computation.
}
\seealso{
  \code{\link{rcomp}}, \code{\link{summary.aplus}}, \code{\link{summary.acomp}}
}
\examples{
data(SimulatedAmounts)
summary(rcomp(sa.lognormals))

}
\keyword{multivariate}
