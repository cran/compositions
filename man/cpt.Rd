\name{cpt}
\alias{cpt}
\alias{cpt.inv}
\title{Centered planar transform}
\description{
  Compute the centered planar  transform of a (dataset of)
  compositions and its inverse.
}
\usage{
          cpt( x )
          cpt.inv( z )
          }
\arguments{
  \item{x}{a composition or a data.matrix of compositions, not necessarily closed}
  \item{z}{the cpt-transform of a composition or a data matrix of
    cpt-transforms of compositions. It is checked that the z sum up to 0.}
}
\value{
  \code{cpt} gives the centered planar transform,
  \code{cpt.inv} gives closed compositions with the given cpt-transforms.
  }
\details{
  The cpt-transform maps a composition in the D-part real-simplex
  isometrically to a D-1 dimensional euclidian vector space, identified with a plane parallel 
  to the simplex but passing through the origin. However the
  transformation is not injective and does not even reach the whole
  plane. Thus resulting covariance matrices are always singular.
  \cr 

  The data can then
  be analysed in this transformed space by all classical multivariate
  analysis tools not relying on a full rank of the covariance matrix. See
  \code{\link{ipt}} and \code{\link{apt}} for alternatives. The
  interpretation of the results is relitevly easy since the relation of each 
  transformed component to the original parts is preserved.\cr
    
  The centered planar transform is given by
  \deqn{ cpt(x)_i := clo(x)_i - \frac1D }{
         cpt(x)_i := clo(x)_i - 1/D}
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{clr}},\code{\link{apt}},\code{\link{ipt}}
}
\references{
}
\examples{
(tmp <- cpt(c(1,2,3)))
cpt.inv(tmp)
cpt.inv(tmp) - clo(c(1,2,3)) # 0
data(Hydrochem)
cdata <- Hydrochem[,6:19]
pairs(cpt(cdata)) 
}
\keyword{multivariate}
