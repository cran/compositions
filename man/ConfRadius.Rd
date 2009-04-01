\name{ConfRadius}
\alias{ConfRadius}
\title{Helper to compute confidence ellipsoids}
\description{
  Computes the quantil of the Mahalanobis distance needed to draw
  confidence ellipsoids.
}
\usage{
ConfRadius(model,prob=1-alpha,alpha)
          }
\arguments{
  \item{model}{A multivariate linear model}
  \item{prob}{The confidence probability}
  \item{alpha}{The alpha error allowed, i.e. the complement of the
    confidence probability}
}
\value{
  a scalar
}
\details{
  Calculates the radius to be used in confidence ellipses for the
  parameters based on the Hottelings \eqn{T^2}{T^2} distribution. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{lm}}, \code{\link{mvar}}, \code{\link{AIC}} 
}
\examples{
data(SimulatedAmounts)
model <- lm(ilr(sa.groups)~sa.groups.area)
vcovAcomp(model)[,,1,1]

}
\keyword{multivariate}
