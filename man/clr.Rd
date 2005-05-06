\name{clr}
\alias{clr}
\alias{clr.inv}
\title{Centered log ratio transform}
\description{
  Compute the centered log ratio transform of a (dataset of)
  composition(s) and its inverse.
}
\usage{
          clr( x )
          clr.inv( z )
          }
\arguments{
  \item{x}{a composition or a data matrix of compositions, not necessarily closed}
  \item{z}{the clr-transform of a composition or a data matrix of
    clr-transforms of compositions, not necessarily centered (i.e. summing up to zero)}
}
\value{
  \code{clr} gives the centered log ratio transform,
  \code{clr.inv} gives closed compositions with the given clr-transforms
  }
\details{
  The clr-transform maps a composition in the D-part Aitchison-simplex
  isometrically to a D-1 dimensonal euclidian vector subspace: consequently, the
  transformation is not injective and only yields vectors which elements
  sum up to 0. Thus resulting covariance matrices are always singular.
  \cr 

  The data can then
  be analysed in this transformation by all classical multivariate
  analysis tools not relying on a full rank of the covariance. See
  \code{\link{ilr}} and \code{\link{alr}} for alternatives. The
  interpretation of the results is relatively easy since the relation between each original
  part and a transformed variable is preserved.\cr
    
  The centered logratio transform is given by
  \deqn{ clr(x) := \left(\ln x_i - \frac1D \sum_{j=1}^D \ln x_j\right)_i }{
    clr(x) := (\bold{ln} \emph{xi} - mean(\bold{ln} x_j)}
  The image of the \code{clr} is given by the vectors with entries
  summing to 0. This hyperplane is also called the clr-plane.

  
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{ilr}},\code{\link{alr}},\code{\link{apt}}
}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.
  

}
\examples{
(tmp <- clr(c(1,2,3)))
clr.inv(tmp)
clr.inv(tmp) - clo(c(1,2,3)) # 0
data(Hydrochem)
cdata <- Hydrochem[,6:19]
pairs(clr(cdata)) 
}
\keyword{multivariate}
