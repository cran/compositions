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
  a generated random dataset of class \code{"acomp"} or \code{"rcomp"},
  drawn from a Dirichlet distribution with the given parameter
  \code{alpha}. The names of \code{alpha} are used to name the parts. 
}
\details{
  The Dirichlet distribution is the result of closing a vector of equally-scaled
  Gamma-distributed variables. It the conjugate prior distribution for a vector
  of probabilities of a multinomial distribution. Thus, it generalizes the beta 
  distribution for more than two parts.
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

  Mateu-Figueras, G.; Pawlowsky-Glahn, V. (2005). The Dirichlet distribution 
  with respect to the Aitchison measure on the simplex, a first approach.
  In: Mateu-Figueras, G. and Barcel\'o-Vidal, C. (Eds.)
  \emph{Proceedings of the 2nd International Workshop on Compositional Data Analysis},
  Universitat de Girona, ISBN 84-8458-222-1, \url{http://ima.udg.es/Activitats/CoDaWork05/}\cr
}
\examples{
tmp <- rDirichlet.acomp(10,alpha=c(A=2,B=0.2,C=0.2))
plot(tmp)
}
\keyword{multivariate}
