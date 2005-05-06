\name{segments}
\alias{segments.rmult}
\alias{segments.acomp}
\alias{segments.rcomp}
\alias{segments.aplus}
\alias{segments.rplus}
\title{Draws straight lines from point to point.}
\description{
  The function draws lines from a points x to a point y in the given geometry.
}
\usage{
          \method{segments}{acomp}(x0,y,...,steps=30)
          \method{segments}{rcomp}(x0,y,...,steps=30)
          \method{segments}{aplus}(x0,y,...,steps=30)
          \method{segments}{rplus}(x0,y,...,steps=30)
          \method{segments}{rmult}(x0,y,...,steps=30)
          }
\arguments{
  \item{x0}{dataset of points of the given type to draw the line from}
  \item{y}{dataset of points of the given type to draw the line to}
  \item{\dots}{further graphical parameters}
  \item{steps}{the number of discretisation points to draw the segments not
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
\seealso{
  \code{\link{plot.acomp}},\code{\link{lines.acomp}} 
}

\examples{
data(SimulatedAmounts)

plot(acomp(sa.lognormals))
segments.acomp(acomp(c(1,2,3)),acomp(c(2,3,1)),col="red")
segments.rcomp(acomp(c(1,2,3)),acomp(c(2,3,1)),col="blue")

plot(aplus(sa.lognormals[,1:2]))
segments.aplus(aplus(c(10,20)),aplus(c(20,10)),col="red")
segments.rplus(rplus(c(10,20)),rplus(c(20,10)),col="blue")

plot(rplus(sa.lognormals[,1:2]))
segments.aplus(aplus(c(10,20)),aplus(c(20,10)),col="red")
segments.rplus(rplus(c(10,20)),rplus(c(20,10)),col="blue")




}
\keyword{multivariate}
