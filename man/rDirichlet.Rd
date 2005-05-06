\name{rDirichlet}
\alias{rDirichlet}
\alias{rDirichlet.acomp}
\alias{rDirichlet.rcomp}
\title{Dirichlet distribution}
\description{
  The Dirichlet distribution on the simplex.
}
\usage{
rDirichlet.acomp(n,alpha)
rDirichlet.rcomp(n,alpha)

          }
\arguments{
  \item{n}{number of datasets to be simulated}
  \item{alpha}{parameters of the Dirichlet distribution}
}
\value{
  a generated random dataset of class \code{"acomp"} or \code{"rcomp"} with
  drawn from a Dirichlet distribution with the given parameter
  \code{alpha}. The names of \code{alpha} are used to name the parts. 
}
\details{
  TO DO!!!
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de},
  Raimon Tolosana-Delgado}
\seealso{
  \code{\link{rnorm.acomp}}
}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr


}
\examples{
tmp <- rDirichlet.acomp(10,alpha=c(A=2,B=0.2,C=0.2))
plot(tmp)
}
\keyword{multivariate}
