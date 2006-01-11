\name{ellipses}
\alias{ellipses}
\alias{ellipses.rmult}
\alias{ellipses.acomp}
\alias{ellipses.rcomp}
\alias{ellipses.aplus}
\alias{ellipses.rplus}
\title{Draw ellipses}
\description{
  Draws ellipses from a mean and a variance into a plot.
}
\usage{
          ellipses(mean,...)
          \method{ellipses}{acomp}(mean,var,r=1,\dots,steps=360)
          \method{ellipses}{rcomp}(mean,var,r=1,\dots,steps=360)
          \method{ellipses}{aplus}(mean,var,r=1,\dots,steps=360)
          \method{ellipses}{rplus}(mean,var,r=1,\dots,steps=360)
          \method{ellipses}{rmult}(mean,var,r=1,\dots,steps=360)
          }
\arguments{
  \item{mean}{a dataset/vector of means or midpoints of the ellipses}
  \item{var}{a variance matrix or a set of variance matrices given by
    \code{var[i,,]}. The principal axis of the variance give the axis of
  the ellipses, whereas the square-root of the eigenvalues times r give the
  half-diameters of the ellipse.}
  \item{r}{a scaling of the half-diameters}
  \item{\dots}{further graphical parameters}
  \item{steps}{the number of discretisation points to draw the ellipses.}
}
\details{
  The ellipse drawn is given by the solutions of
  \deqn{(x-mean)^tvar^{-1}(x-mean)=r^2}
  in the respective geometry of the parameter space. Note that
  these ellipses can be added to panel plots (by means of 
  orthogonal projections in the corresponding geometry).
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{plot.acomp}}, 
}

\examples{
data(SimulatedAmounts)

plot(acomp(sa.lognormals))
tt<-acomp(sa.lognormals); ellipses(mean(tt),var(tt),r=2,col="red")
tt<-rcomp(sa.lognormals); ellipses(mean(tt),var(tt),r=2,col="blue")

plot(aplus(sa.lognormals[,1:2]))
tt<-aplus(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="red")
tt<-rplus(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="blue")

plot(rplus(sa.lognormals[,1:2]))
tt<-aplus(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="red")
tt<-rplus(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="blue")
tt<-rmult(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="green")

}
\keyword{multivariate}
