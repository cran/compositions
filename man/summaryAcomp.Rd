\name{summary.acomp}
\alias{summary.acomp}
\title{Summarizing a compositional dataset in terms of ratios}
\description{
  Summaries in terms of compositions are quite different from classical
  ones. Instead of analysing each variable individually, we must
  analyse each pairwise ratio in a log geometry.
}
\usage{
          \method{summary}{acomp}( object, \dots )
          }
\arguments{
  \item{object}{a data.matrix of compositions, not necessarily closed}
  \item{\dots}{not used, only here for generics}
}
\value{
  The result is an object of type \code{"summary.acomp"}
  \item{mean}{The \code{\link{mean.acomp}} composition}
  \item{mean.ratio}{A matrix containing the geometric mean of the
    pairwise ratios}
  \item{variation}{The variation matrix of the dataset (\code{{\link{variation.acomp}}})}
  \item{expsd}{A matrix containing the  one-sigma factor for
    each ratio, computed as \code{exp(sqrt(variation.acomp(W)))}. To
    obtain two-sigma-factor it needs to be squared. To obtain the
    reverse bound we compute 1/expsd}
  \item{min}{A matrix containing the minimum of each of the pairwise ratios}
  \item{q1}{A matrix containing the 1-Quartile of each of the pairwise ratios}
  \item{median}{A matrix containing the median of each of the pairwise ratios}
  \item{q1}{A matrix containing the 3-Quartile of each of the pairwise ratios}
  \item{max}{A matrix containing the maximum of each of the pairwise
    ratios}

}
\details{
  It is quite difficult to summarize a composition in a consistent and
  interpretable way. We tried to provide such a summary here.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{acomp}}
}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.
}
\examples{
data(SimulatedAmounts)
summary(acomp(sa.lognormals))

}
\keyword{multivariate}
