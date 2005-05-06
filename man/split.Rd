\name{split}
\alias{split.acomp}
\alias{split.aplus}
\alias{split.rcomp}
\alias{split.rplus}
\alias{split.rmult}
\title{Spliting datasets in groups given by factors}
\description{
  Each of the given scale levels has an associated norm, which is
  computed for each element by these functions.
}
\usage{
\method{split}{acomp}(x,f)
\method{split}{rcomp}(x,f)
\method{split}{aplus}(x,f)
\method{split}{rplus}(x,f)
\method{split}{rmult}(x,f)
          }
\arguments{
  \item{x}{a dataset or a single vector of some type}
  \item{f}{a factor that defines the grouping or a list of factors }
}
\value{
  a list of objects of the same type as x.
}
\examples{
  data(SimulatedAmounts)
  split(acomp(sa.groups),sa.groups.area)
}
\seealso{
  \code{\link[base]{split}}
}
\keyword{}
