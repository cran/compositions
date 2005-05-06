\name{gsi.drop}
\alias{gsi.drop}
\title{Internal functions: A conditional drop}
\description{
  drop, if drop is needed.
}
\usage{
gsi.drop(X,drop)
}
\arguments{
  \item{X}{ an array needing dimensions dropped}
  \item{drop}{ a logical whether to drop dimensions}
}
\value{
  \code{X} or \code{drop(X)}
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\details{

}

\note{Do not use gsi.* functions directly since they are internal
  functions of the package}

\seealso{
  \code{\link[base]{drop}}, \code{\link{gsi}}
}
\examples{

}
\keyword{}
