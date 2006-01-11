\name{ult}
\alias{ult}
\alias{ult.inv}
\alias{Kappa}
\title{Uncentered log transform}
\description{
  Compute the uncentered log ratio transform of a (dataset of)
  composition(s) and its inverse.
}
\usage{
          ult( x ,...)
          ult.inv( z ,...)
          Kappa( x ,...)
          }
\arguments{
  \item{x}{a composition or a data matrix of compositions, not necessarily closed}
  \item{z}{the ult-transform of a composition or
    clr-transforms of compositions (or a data matrix), not necessarily centered}
  \item{\dots}{for generic use only}

}
\value{
  \code{ult} gives the uncentered log transform,\cr
  \code{ult.inv} gives closed compositions with the given
  ult/clr-transforms\cr
  \code{Kappa} gives the difference between the clr and the ult
  transforms. It is quite linked to information measures.
  }
\details{
  The ult-transform is simply the elementwise log of the closed
  composition. The ult has some important properties in the scope 
  of Information Theory of probability vectors (but might be mostly 
  misleading for exploratory analysis of compositions).

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{ilr}},\code{\link{alr}},\code{\link{apt}}
}
\examples{
(tmp <- ult(c(1,2,3)))
ult.inv(tmp)
ult.inv(tmp) - clo(c(1,2,3)) # 0
data(Hydrochem)
cdata <- Hydrochem[,6:19]
pairs(ult(cdata))
Kappa(c(1,2,3))
}
\keyword{multivariate}
