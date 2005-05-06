\name{gsi.simshape}
\alias{gsi.simshape}
\title{Internal function: Reshape an object to the shape type of another}
\description{
Reshape an object to the shape type of another
}
\usage{
          gsi.simshape(x,oldx)
          }
\arguments{
  \item{x}{the data to be returned}
  \item{oldx}{a data of the intended shape}
}
\value{
  x is shaped as oldx
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\details{
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
}
\references{
}
\examples{
gsi.simshape(matrix(1:4,nrow=1),1:3)
}
\keyword{}
