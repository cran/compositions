\name{coloredBiplot}
\alias{coloredBiplot}
\alias{coloredBiplot.default}
\alias{coloredBiplot.princomp}
\alias{coloredBiplot.prcomp}
\title{A biplot providing somewhat easier access to details of the plot.}
\description{
  This function generates a simple biplot out of various source and
  allows to color the objects individually. (Warning: This preliminary
  version is likely to be
  changed in next edition)
}
\usage{

\method{coloredBiplot}{default}(x, y, var.axes = TRUE, col, 
         cex = rep(par("cex"), 2), xlabs = NULL, ylabs = NULL, expand=1, 
         xlim = NULL, ylim = NULL, arrow.len = 0.1, main = NULL, sub = NULL, 
         xlab = NULL, ylab = NULL, xlabs.col = NULL, xlabs.pc=NULL, ...)
\method{coloredBiplot}{princomp}(x, choices = 1:2, scale = 1, 
         pc.biplot=FALSE, ...)
\method{coloredBiplot}{prcomp}(x, choices = 1:2, scale = 1, 
         pc.biplot=FALSE, ...)
         
}
\arguments{
  \item{x}{a representation of the the co-information to be
    plotted, given by a result of princomp or prcomp; or the first set of
    coordinates to be plotted}
  \item{y}{optional, the second set of coordinates to be potted}
  \item{var.axes}{If 'TRUE# the second set of points have arrows representing
    them as (unscaled) axes.}
  \item{col}{One color (to be used for the y set) or a vector of two colors 
    (to be used for x and y sets respectively, if xlabs.col is NULL)}
  \item{cex}{the usual cex parameter for plotting; can be a length-2 vector to
    format differently x and y labels/symbols.}
  \item{xlabs}{The names to write for the points of the first set}
  \item{ylabs}{The names to write for the points of the second set}
  \item{expand}{An expansion factor to apply when plotting the second set of 
    points relative to the first. This can be used to tweak the scaling of the 
    two sets to a physically comparable scale.}
  \item{xlim}{horizontal axis limits}  
  \item{ylim}{vertical axis limits}  
  \item{arrow.len}{The length of the arrow heads on the axes plotted in 
  'var.axes' is true. The arrow head can be suppressed by 'arrow.len=0'.}
  \item{main}{main title}  
  \item{sub}{subtitle}  
  \item{xlab}{horizontal axis title}  
  \item{ylab}{vertical axis title}  
  \item{xlabs.col}{The color(s) to draw the points of the first set}
  \item{xlabs.pc}{the plotting character(s) for the first set} 
  \item{scale}{the way to distribute the singular values on the 
    right or left singular vectors for princomp and prcomp objects 
    (see \code{\link{biplot}})}
  \item{choices}{the components to be plotted (see \code{\link{biplot}})}
  \item{pc.biplot}{should be scaled by \code{sqrt(nrow(X))}? 
    (see \code{\link{biplot}}) }
  \item{\dots}{further parameters for plot}
}
\value{
The function is called only for the side effect of plotting. It is a modification 
of the standard R routine 'biplot'.
}
\details{
  The functions is provided for convenience. 
}
\author{Raimon Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{biplot}}, \code{\link{plot.acomp}}
}
\examples{
data(SimulatedAmounts)
coloredBiplot(x=princomp(acomp(sa.outliers5)),pc.biplot=FALSE,
          xlabs.pc=c(1,2,3), xlabs.col=2:4, col="black")
}
\keyword{multivariate}
