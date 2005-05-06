\name{runif}
\alias{runif.acomp}
\alias{runif.rcomp}
\title{The uniform distribution on the simplex}
\description{
  Generates random compositions with a uniform distribution on the
  (rcomp) simplex. 
}
\usage{
runif.acomp(n,D)
runif.rcomp(n,D)
          }
\arguments{
  \item{n}{number of datasets to be simulated}
  \item{D}{number of parts}
}
\value{
  a generated random dataset of class \code{"acomp"} or \code{"rcomp"} with
  drawn from a uniform distribution on the simplex of D parts. }
\seealso{
  \code{\link{rDirichlet.acomp}}
}
\note{The only difference between both routines is the type of
the returned dataset.}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr
}
\examples{
plot(runif.acomp(10,3))
plot(runif.rcomp(10,3))
}
\keyword{multivariate}
