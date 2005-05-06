\name{variation}
\alias{variation}
\alias{variation.default}
\alias{variation.acomp}
\alias{variation.rcomp}
\alias{variation.aplus}
\alias{variation.rplus}
\alias{variation.rmult}
\title{Variation matrices of amounts and compositions}
\description{
  Compute the variation  matrix in the various approaches of compositional
 and amount  data analysis. Pay attention that this is not computing the variance or
 covariance matrix!
}
\usage{
          variation(x,\dots)
          \method{variation}{acomp}(x, \dots)
          \method{variation}{rcomp}(x, \dots)
          \method{variation}{aplus}(x, \dots)
          \method{variation}{rplus}(x, \dots)
          \method{variation}{rmult}(x, \dots)
          }
\arguments{
  \item{x}{a dataset, eventually of amounts or compositions}
  \item{\dots}{currently unused}
}
\value{
  The variation matrix of x.
}
\details{
   The variation matrix was defined in the acomp context for 
   analysis of compositions as the matrix of variances of all 
   possible log-ratios among components (Aitchison, 1986). The 
   generalization to rcomp objects is simply to reproduce the 
   variance of all possible differences between components. The 
   amount and rmult objects should not be treated with variation 
   matrices, because this implies always the existence of a closure.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{cdt}}, \code{\link{clrvar2ilr}}, \code{\link{clo}},
\code{\link{mean.acomp}}, \code{\link{acomp}}, \code{\link{rcomp}},
\code{\link{aplus}}, \code{\link{rplus}} }

\examples{
data(SimulatedAmounts)
mean.col(sa.lognormals)
variation(acomp(sa.lognormals))
variation(rcomp(sa.lognormals))
variation(aplus(sa.lognormals))
variation(rplus(sa.lognormals))
variation(rmult(sa.lognormals))

}
\keyword{multivariate}
