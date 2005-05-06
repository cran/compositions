\name{boxplot}
\alias{boxplot.acomp}
\alias{boxplot.rcomp}
\alias{boxplot.rplus}
\alias{boxplot.aplus}
\alias{vp.boxplot}
\alias{vp.logboxplot}
\title{Displaying compositions and amounts by boxplots}
\description{
  For the different interpretations of amount data a different type of
  boxplot is feasible. Thus different boxplots are drawn.
}
\usage{
\method{boxplot}{acomp}(x,fak=NULL,\dots,
                         xlim=x.lim,ylim=c(minq,maxq),log=TRUE,
                         panel=vp.logboxplot,dots=!boxes,boxes=TRUE)
\method{boxplot}{rcomp}(x,fak=NULL,\dots,
                         xlim=x.lim,ylim=c(0,1),log=FALSE,
                         panel=vp.boxplot,dots=!boxes,boxes=TRUE)
\method{boxplot}{aplus}(x,fak=NULL,\dots,log=TRUE)
\method{boxplot}{rplus}(x,fak=NULL,\dots,ylim=c(0,max(x)),log=FALSE)
vp.boxplot(x,y,...,dots=FALSE,boxes=TRUE,xlim,ylim,log,notch=FALSE) 
vp.logboxplot(x,y,...,dots=FALSE,boxes=TRUE,xlim,ylim,log,notch=FALSE) 
          }
\arguments{
  \item{x}{a dataset}
  \item{fak}{a factor to split the dataset, not yet implemented in aplus
  and rplus}
  \item{xlim}{x-limits of the plot}
  \item{ylim}{y-limits of the plot}
  \item{log}{logical indicating whether ploting should be done on log
    scale}
  \item{panel}{the panel function to be used or a list of multiple panel
    functions}
  \item{\dots}{further graphical parameters}
  \item{dots}{a logical indicating whether the points should be drawn}  
  \item{boxes}{a logical indicating the boxes should be drawn}
  \item{y}{used by pairs}
  \item{notch}{should the boxes be notched}
}
\details{
  \code{vp.boxplot} and \code{vp.logboxplot} are only used as panel functions.
}
\seealso{
  \code{\link{plot.acomp}}, \code{\link{qqnorm.acomp}} 
}

\examples{
data(SimulatedAmounts)
boxplot(acomp(sa.lognormals))
boxplot(rcomp(sa.lognormals))
boxplot(aplus(sa.lognormals))
boxplot(rplus(sa.lognormals))
}
\keyword{multivariate}
