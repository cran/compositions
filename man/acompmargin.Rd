\name{acompmargin}
\alias{acompmargin}
\title{Marginal compositions in Aitchison Compositions}
\description{
  Compute marginal compositions of selected parts, by computing the rest as the
  geometric mean of the non-selected parts.
}
\usage{
          acompmargin(X,d=c(1,2),name="*",pos=length(d)+1)
          }
\arguments{
  \item{X}{composition or dataset of compositions}
  \item{d}{vector containing the indices xor names of the columns selected}
  \item{name}{The new name of the amalgamation column}
  \item{pos}{The position where the new amalgamation column should be
    stored. This defaults to the last column.}
}
\value{
  A closed compositions with class \code{"acomp"} containing the
  variables given by \code{d} and the the amalgamation column.  }
\details{ The amalgamation column is simply computed by taking the
  geometric mean of the non-selected components. This is
consistent with the \code{\link{acomp}} approach and gives clear ternary
diagrams. However, this geometric mean is difficult to interpret.
}
\seealso{ \code{\link{rcompmargin}}, \code{\link{acomp}}
}
\references{
  Vera Pawlowsky-Glahn (2003) personal communication. Seminaris 
  d'Estadística, Departament d'Informàtica i Matemàtica Aplicada,
  Universitat de Girona.
}
\examples{
data(SimulatedAmounts)
plot.acomp(sa.lognormals5,margin="acomp")
plot.acomp(acompmargin(sa.lognormals5,c("Pb","Zn")))
plot.acomp(acompmargin(sa.lognormals5,c(1,2)))
}
\keyword{multivariate}
