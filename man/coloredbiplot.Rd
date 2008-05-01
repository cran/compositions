\name{coloredBiplot}
\alias{coloredBiplot}
\title{A biplot providing somewhat easier access to details of the plot.}
\description{
  This function generates a simple biplot out of various source and
  allows to color the objects individually. (Warning: This preliminary
  version is likely to be
  changed in next edition)
}
\usage{
coloredBiplot(xrf,xcol="black",ycol="red",scale=1,choice=c(1,2),pc.biplot=FALSE,pch=19,ypch=2,cex=1,xarrows=FALSE,yarrows=TRUE,...)

}
\arguments{
  \item{xrf}{some representation of the the coinformation to be
    plotted given by a result of princomp, prcomp or svd}
  \item{xcol}{The color(s) to draw the points of the first set}
  \item{ycol}{The color(s) to draw points or arrows of second set}
  \item{scale}{exponential scaling ???}
  \item{choice}{the components to be plotted}
  \item{pc.biplot}{should be scaled by \code{sqrt(nrow(X))} ???}
  \item{pch}{the plotchar for the first set} 
  \item{ypch}{the plotchar for the y component}
  \item{cex}{the usual cex parameter for plotting}
  \item{xarrows}{wether or not to plot the arrows for the first dataset}  
  \item{yarrows}{wether or not to plot the arrows for the second dataset}  
  \item{\dots}{further parameters for plot}
}
\value{
The function is called only for the sideeffect of plotting. 
}
\details{
  The functions is provided for conveniance and is to be considered
experimental. 
}
\note{This function is experimental a more elaborated incorporating
generics and more detailed controls will be provided in the next version
or a seperate package.}
\author{Raimon Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{biplot}}, \code{\link{plot.acomp}}
}

\examples{
data(SimulatedAmounts)
coloredBiplot(xrf=svd(var(sa.outliers5)),pc.biplot=FALSE,pch=c(1,2,3),xcol=1:3,ycol=c("blue","cyan","magenta"))

}
\keyword{multivariate}
