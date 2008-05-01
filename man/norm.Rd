\name{norm}
%\alias{norm}
\alias{norm.default}
\alias{norm.acomp}
\alias{norm.aplus}
\alias{norm.rcomp}
\alias{norm.rplus}
\alias{norm.rmult}
\title{Vector space norm}
\description{
  Each of the considered space structures has an associated norm, which is
  computed for each element by these functions.
}
\usage{
\method{norm}{default}(X,\dots)
\method{norm}{acomp}(X,\dots)
\method{norm}{rcomp}(X,\dots)
\method{norm}{aplus}(X,\dots)
\method{norm}{rplus}(X,\dots)
\method{norm}{rmult}(X,\dots)

          }
\arguments{
  \item{X}{a dataset or a single vector of some type}
  \item{\dots}{currently not used, intended to select a 
              different norm rule in the future}
}
\value{
  The norms of the given vectors.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{normalize}}
}
\examples{
data(SimulatedAmounts)
tmp <- acomp(sa.lognormals)
mvar(tmp)
sum(norm( tmp - mean(tmp) )^2)/(nrow(tmp)-1)

}
\keyword{multivariate}
