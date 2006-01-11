\name{barplot.acomp}
\alias{barplot.acomp}
\alias{barplot.rcomp}
\alias{barplot.aplus}
\alias{barplot.rplus}
\title{Bar charts of amounts}
\description{
  Compositions and amounts dispalyed as bar plots. 
}
\usage{
\method{barplot}{acomp}(height,...,legend.text=TRUE,beside=FALSE,total=1)
\method{barplot}{rcomp}(height,...,legend.text=TRUE,beside=FALSE,total=1)
\method{barplot}{aplus}(height,...,legend.text=TRUE,beside=TRUE)
\method{barplot}{rplus}(height,...,legend.text=TRUE,beside=TRUE)
          }
\arguments{
  \item{height}{an acomp, rcomp, aplus, or rplus object giving amounts
    to be displayed}
  \item{\dots}{further graphical parameters as in
    \code{\link[graphics]{barplot}}}
  \item{legend.text}{same as legend.text in \code{\link[graphics]{barplot}}}
  \item{beside}{same as beside in \code{\link[graphics]{barplot}}}
  \item{total}{The total to be used in displaying the composition,
    typically 1, 100 or the number of parts}
}
\value{A numeric vector (or matrix, when \code{beside = TRUE}) giving
  the coordinates of all the bar midpoints drawn, as in
  \code{\link[graphics]{barchart}} }
\details{
  The functions are essentially light-weighted wrappers for
  \code{\link[graphics]{barplot}}, just adding an adequate default
  behavior for each of the scales.  


}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{acomp}}, \code{\link{rcomp}},   \code{\link{rplus}}
  \code{\link{aplus}}, \code{\link{plot.acomp}},
  \code{\link{boxplot.acomp}}


}
\examples{
data(SimulatedAmounts)
barplot(mean(acomp(sa.lognormals[1:10,])))
barplot(mean(rcomp(sa.lognormals[1:10,])))
barplot(mean(aplus(sa.lognormals[1:10,])))
barplot(mean(rplus(sa.lognormals[1:10,])))

barplot(acomp(sa.lognormals[1:10,]))
barplot(rcomp(sa.lognormals[1:10,]))
barplot(aplus(sa.lognormals[1:10,]))
barplot(rplus(sa.lognormals[1:10,]))

}
\keyword{multivariate}
