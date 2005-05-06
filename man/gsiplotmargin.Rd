\name{gsiinternal}
\alias{gsi.plotmargin}
\title{Internal functions of the compositions package}
\description{
  Internal function to compute 2D marginal compositions for plots  
}
\usage{
gsi.plotmargin(X,d,margin)
}
\arguments{
  \item{X}{A multivariate compositional dataset}
  \item{d}{a numeric or character vector of two elements specifying the
    margin}
  \item{margin}{a character specifying the type of margin to be
    choosen. Possible values are "acomp", "rcomp" or a column name from
    the dataset.}
}
\value{
  a composition of three elements.
}
\details{

}

\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}
\seealso{
  \code{\link{gsi}}
}
\examples{
}
\keyword{}
