\name{lines}
\alias{lines.rmult}
\alias{lines.acomp}
\alias{lines.rcomp}
\alias{lines.aplus}
\alias{lines.rplus}
\title{Draws connected lines from point to point.}
\description{
  Functions taking coordinates given in various ways and
     joining the corresponding points with line segments.
}
\usage{
          \method{lines}{acomp}(x,...,steps=30)
          \method{lines}{rcomp}(x,...,steps=30)
          \method{lines}{aplus}(x,...,steps=30)
          \method{lines}{rplus}(x,...,steps=30)
          \method{lines}{rmult}(x,...,steps=30)
          }
\arguments{
  \item{x}{a dataset of the given type}
  \item{\dots}{further graphical parameters}
  \item{steps}{the numer of discretisation points to draw the segments not
    straight on the monitor.}
}
\details{
The functions add lines to the graphics generated with the corresponding
plot functions.
\cr

Adding to multipaneled plots, redraws the plot completely and is only
possible, when the plot has been created with the plotting routines from
this library.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}
\seealso{
  \code{\link{plot.acomp}}, \code{\link{straight}} 
}

\examples{
data(SimulatedAmounts)

plot(acomp(sa.lognormals))
lines(acomp(sa.lognormals),col="red")
lines(rcomp(sa.lognormals),col="blue")

plot(aplus(sa.lognormals[,1:2]))
lines(aplus(sa.lognormals[,1:2]),col="red")
lines(rplus(sa.lognormals)[,1:2],col="blue")

plot(rplus(sa.lognormals[,1:2]))
tt<-aplus(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="red")
tt<-rplus(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="blue")
tt<-rmult(sa.lognormals[,1:2]); ellipses(mean(tt),var(tt),r=2,col="green")

}
\keyword{multivariate}
