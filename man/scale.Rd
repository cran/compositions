\name{scale}
\alias{scale.acomp}
\alias{scale.aplus}
\alias{scale.rcomp}
\alias{scale.rplus}
\alias{scale.rmult}
\title{Normalizing datasets by centering and scaling}
\description{
  The dataset is standardized by optional scaling and centering.
}
\usage{
\method{scale}{acomp}(x,center=TRUE, scale=TRUE)
\method{scale}{rcomp}(x,center=TRUE, scale=TRUE)
\method{scale}{aplus}(x,center=TRUE, scale=TRUE)
\method{scale}{rplus}(x,center=TRUE, scale=TRUE)
\method{scale}{rmult}(x,center=TRUE, scale=TRUE)
          }
\arguments{
  \item{x}{a dataset or a single vector of some type}
  \item{center}{logical value}
  \item{scale}{logical value}
}
\value{
  a vector or data matrix, as x and with the same class, but acordingly transformed.
}
\details{
  scaling is defined in various ways for the different data types. It is
  always performed as an operation in the enclosing vector space. Not in
  all cases an independent scaling of the different coordinates is
  appropriate. This is only done for rplus and rmult.
  }
\examples{
  data(SimulatedAmounts)
  plot(scale(acomp(sa.groups)))
  \dontrun{
plot(scale(rcomp(sa.groups)))
}
  plot(scale(aplus(sa.groups)))
  \dontrun{
plot(scale(rplus(sa.groups)))
}
  plot(scale(rmult(sa.groups)))

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{split}{base}}
}
\keyword{multivariate}
