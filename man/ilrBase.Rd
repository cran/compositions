\name{ilrBase}
\alias{ilrBase}
\alias{gsi.ilrBase}
\alias{ilrBaseList}
\title{The canonical basis in the clr plane used for ilr and ipt transforms.}
\description{
  Compute the basis of a clr-plane, to use with isometric log-ratio or planar transform of a (dataset of)
  compositions.
}
\usage{
          ilrBase( x=NULL , z=NULL , D = NULL )
          gsi.ilrBase(D)
          }
\arguments{
  \item{x}{ optional dataset or vector of compositions}
  \item{z}{ optional dataset or vector containing ilr or ipt coordinates} 
  \item{D}{ number of parts of the simplex}
}
\value{
  Both methods give a matrix containing by columns the basis elements for the
  canonical basis of the clr-plane used for the ilr and ipt transform. 
}
\details{
  ilrBase is a wrapper catching the answers of gsi.ilrBase and is to be
  used as the more convenient function. Only one of the arguments is needs 
  to determine the dimension of the simplex.
}
\seealso{
  \code{\link{clr}},\code{\link{ilr}},\code{\link{ipt}},
  \url{http://ima.udg.es/Activitats/CoDaWork03}
}
\references{
  Egozcue J.J., V. Pawlowsky-Glahn, G. Mateu-Figueras and
  C. Barcel'o-Vidal (2003) Isometric logratio transformations for
  compositional data analysis. \emph{Mathematical Geology}, \bold{35}(3)
  279-300\cr
  }
\examples{
ilr(c(1,2,3))
ilrBase(D=2)
ilrBase(c(1,2,3))
ilrBase(z= ilr(c(1,2,3)) )
round(ilrBase(D=7),digits= 3)
}
\keyword{multivariate}
