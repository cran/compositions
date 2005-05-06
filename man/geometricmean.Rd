\name{geometricmean}
\alias{geometricmean}
\alias{geometricmean.row}
\alias{geometricmean.col}

\title{The geometric mean}
\description{
  Computes the geometric mean.
}
\usage{
          geometricmean(x,\dots)
          geometricmean.row(x,\dots)
          geometricmean.col(x,\dots)
          
          }
\arguments{
  \item{x}{a numeric vector or matrix of data }
  \item{\dots}{further arguments to compute the mean }
}
\value{
  The geometric means of x as a whole (geometricmean), its rows
  (geometricmean.row) or its columns (geometricmean.col).
}
\details{
  The geometric mean is defined as:
  \deqn{geometricmean(x) := \left( \prod_{i=1}^n x_i\right)^{1/n}}{
    geometricmean(x) := \code{prod(x)^(1/length(x))}
  }
  The geometric mean is actually computed by
  \code{exp(mean(log(c(unclass(x))),...))}.  }

\seealso{
  \code{\link{mean.rplus}} }

\examples{
geometricmean(1:10)

}
\keyword{univar}
