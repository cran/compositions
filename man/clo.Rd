\name{clo}
\alias{clo}
\title{Closure of a composition}
\description{
  Closes compositions to sum up to one (or an optional total), by dividing each part by the sum.
}
\usage{
          clo( X, parts=1:NCOL(oneOrDataset(X)),total=1)
          }
\arguments{
  \item{X}{composition or dataset of compositions}
  \item{parts}{vector containing the indices xor names of the columns to be
    used}
  \item{total}{the total amount to which the compositions should be closed; either 
    a single number, or a numeric vector of length
    \code{gsi.getN(X)} specifying a different total for each compositional vector in the dataset.}
}
\value{
  a composition or a datamatrix of compositions without compositional
  class. The individual compositions are forced to sum to 1 (or to 
  the optionally-specified total). The result
  should have the same shape as the input (vector, row, matrix). 
}
\details{
  The closure operation is given by
  \deqn{ clo(x) := \left(x_i / \sum_{j=1}^D  x_j\right) }{
    clo(x) := ( \emph{x_i} / sum(x_j))}
  \code{clo} makes a composition without assigning one of the
  compositional classes \code{\link{acomp}} or \code{\link{rcomp}}. 
  Note that after computing the closed-to-one version, obtaining a 
  version closed to any other value is done by simple multiplication.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}

\seealso{
  \code{\link{clr}},\code{\link{acomp}},\code{\link{rcomp}}
}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.
  

}
\note{\code{clo} can be used to unclass compositions.}
\examples{
(tmp <- clo(c(1,2,3)))
clo(tmp,total=100)
data(Hydrochem)
cdata <- Hydrochem[,6:19]
plot( clo(Hydrochem,8:9) ) # Giving points on a line 

}
\keyword{multivariate}
