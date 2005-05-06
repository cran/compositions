\name{ipt}
\alias{ipt}
\alias{ipt.inv}
\alias{ucipt.inv}
\title{Isometric planar transform}
\description{
  Compute the isometric planar  transform of a (dataset of)
  composition(s) and its inverse.
}
\usage{
          ipt( x , V = ilrBase(x) )
          ipt.inv( z , V = ilrBase(z=z) )
          ucipt.inv( z , V = ilrBase(z=z) )
          }
\arguments{
  \item{x}{a composition or a data matrix of compositions, not necessarily closed}
  \item{z}{the ipt-transform of a composition or a data matrix of
    ipt-transforms of compositions}
  \item{V}{a matrix with columns giving the chosen basis of the clr-plane}

}
\value{
  \code{ipt} gives the centered planar transform,
  \code{ipt.inv} gives closed compositions with with the given ipt-transforms,
  \code{ucipt.inv} unconstrained ipt.inv does the same as ipt.inv but
  sets illegal values to NA rather then giving an error. This is a
  workaround to allow procedures not honoring the constraints of the
  space.  
  }
\details{
  The ipt-transform maps a composition in the D-part real-simplex
  isometrically to a D-1 dimensonal euclidian vector. Although the
  transformation does not reach the whole \eqn{R^{D-1}}, resulting covariance
  matrices are typically of full rank.
  \cr 
  The data can then
  be analysed in this transformation by all classical multivariate
  analysis tools. However, interpretation of results may be
  difficult, since the
  transform does not keep the variable names, given that there is no
  one-to-one relation between the original parts and each transformed variables. See
  \code{\link{cpt}} and \code{\link{apt}} for alternatives. \cr
    
  The isometric planar transform is given by
    \deqn{ ipt(x) := V^t cpt(x)  }
  with \code{\link{cpt}}(x) the centred planar transform and
  \eqn{V\in R^{d \times (d-1)}}{V} a matrix which columns form an orthonormal 
  basis of the clr-plane. A default matrix \eqn{V} is given by
  \code{ilrBase(\var{D})}

}
\seealso{
  \code{\link{ilr}},\code{\link{ilrBase}}, \code{\link{cpt}}
}
\references{
}
\examples{
(tmp <- ipt(c(1,2,3)))
ipt.inv(tmp)
ipt.inv(tmp) - clo(c(1,2,3)) # 0
data(Hydrochem)
cdata <- Hydrochem[,6:19]
pairs(ipt(cdata)) 
}
\keyword{multivariate}
