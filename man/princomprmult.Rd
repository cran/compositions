\name{princomp.rmult}
\alias{princomp.rmult}
\title{Principal component analysis for real data}
\description{
  Performs a principal component analysis for datasets of type rmult.
}
\usage{
\method{princomp}{rmult}(x,\dots)
}
\arguments{
  \item{x}{a rmult-dataset}
  \item{\dots}{Further arguments to call \code{\link{princomp.default}}}
}
\value{
  An object of type \code{princomp} with the following fields
  \item{sdev}{the standard deviation of the principal components.}
  \item{loadings}{the matrix of variable loadings (i.e., a matrix whose
    columns contain the eigenvectors). This is of class
    \code{"loadings"}.}
  \item{center}{the mean that was substracted from the data set}
  \item{scale}{the scaling applied to each variable}
  \item{n.obs}{number of observations}
  \item{scores}{if \code{scores = TRUE}, the scores of the supplied data
    on the principal components. Scores are coordinates in a basis given by the 
    principal components.}
  \item{call}{the matched call}
  \item{na.action}{Not clearly understood}
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\details{
  The function just does \code{princomp(unclass(x),\dots,scale=scale)}
  and is only here for convenience. 
}
\seealso{
  \code{\link{princomp.rplus}}
}
\examples{
data(SimulatedAmounts)
pc <- princomp(rmult(sa.lognormals5))
pc
summary(pc)
plot(pc) 
screeplot(pc)
screeplot(pc,type="l")
biplot(pc)
biplot(pc,choice=c(1,3))
loadings(pc)
plot(loadings(pc))
pc$sdev^2
cov(predict(pc,sa.lognormals5))
}
\keyword{multivariate}
