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
\method{split}{acomp}(x,f,drop=FALSE,\dots)
\method{split}{rcomp}(x,f,drop=FALSE,\dots)
\method{split}{aplus}(x,f,drop=FALSE,\dots)
\method{split}{rplus}(x,f,drop=FALSE,\dots)
\method{split}{rmult}(x,f,drop=FALSE,\dots)
          }
\arguments{
  \item{x}{a dataset or a single vector of some type}
  \item{f}{a factor that defines the grouping or a list of factors }
  \item{drop}{drop=FALSE also gives (empty) datsets for empty
    categories}
  \item{\dots}{Further arguments passed to split.default.
    Currently probably without any use.}
}
\value{
  a list of objects of the same type as x.
}
\examples{
  data(SimulatedAmounts)
  split(acomp(sa.groups),sa.groups.area)
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link[base]{split}}
}
\keyword{}
