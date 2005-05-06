\name{ilt}
\alias{ilt}
\alias{ilt.inv}
\title{Isometric log transform}
\description{
  Compute the isometric log  transform of a vector (dataset) of
  amounts and its inverse.
}
\usage{
          ilt( x )
          ilt.inv( z  )
          }
\arguments{
  \item{x}{a vector or data matrix of amounts}
  \item{z}{the ilt-transform of a vector or data matrix of
    ilt-transforms of amounts}
}
\value{
  \code{ilt} gives the isometric log transform, i.e. simply the log of
  the argument,
  \code{ilt.inv} gives amounts with the given ilt, i.e. simple the exp
  of the argument
  }
\details{
  The ilt-transform maps  D amounts (considered in log geometry)
  isometrically to a D dimensional euclidean vector. The \code{ilt} is
  part of the \code{\link{aplus}} framework.
  \cr 
  The data can then be analysed in this transformation by all classical
  multivariate analysis tools. The interpretation of the results is easy
  since the relation to the original
  variables is preserved. \cr
    
  The isometric log transform is given by
    \deqn{ ilt(x)_i := \ln x_i }

}
\seealso{
  \code{\link{ilr}}, \code{\link{iit}}, \code{\link{aplus}}
}
\references{
}
\examples{
(tmp <- ilt(c(1,2,3)))
ilt.inv(tmp)
ilt.inv(tmp) - c(1,2,3) # 0
data(Hydrochem)
cdata <- Hydrochem[,6:19]
pairs(ilt(cdata)) 
}
\keyword{multivariate}
