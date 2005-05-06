\name{plot.aplus}
\alias{plot.aplus}
\alias{plot.rplus}
\alias{plot.rmult}
\title{Displaying amounts in scatterplots}
\description{
}
\usage{
\method{plot}{aplus}(x,...,labels=colnames(X),cn=colnames(X),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),add=FALSE,logscale=TRUE,col=par("col"))
\method{plot}{rplus}(x,...,labels=colnames(X),cn=colnames(X),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),add=FALSE,logscale=FALSE,xlim=apply(X,2,function(x) c(0,max(x))),ylim=xlim,col=par("col"))
\method{plot}{rmult}(x,...,labels=colnames(X),cn=colnames(X),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),add=FALSE,logscale=FALSE,col=par("col"))
          }
\arguments{
  \item{x}{a dataset of an amount class}
  \item{\dots}{further graphical parameters passed (see
    \code{\link[graphics]{par}})}
  \item{add}{a logical indicating whether the information should just
    be added to an existing plot. In case of false a new plot is
    created.}
  \item{col}{The color to plot the data.}
  \item{labels}{The names of the parts}
  \item{cn}{The names of the parts to be used in a single
    panel. Internal use only.}
  \item{aspanel}{Logical indicating that only a single panel should be
    drawn and not the whole plot. Internal use only.}
  \item{id}{A logical. If true one can identify the points like with the
    \code{\link[graphics]{identify}} command. }
  \item{idlabs}{A character vector providing the labels to be used with
    the identification, when \code{id=TRUE}}
  \item{idcol}{color of the idlabs-labels}
  \item{center}{a logical indicating whether a the data should be
    centered prior to the plot. Centering is done in the choosen
    philosophy. See \code{\link{scale}}}
  \item{scale}{a logical indicating whether a the data should be
    scaled prior to the plot. Scaling is done in the choosen
    philosophy. See \code{\link{scale}}}
  \item{pca}{A logical indicating whether the first principle component
    should be displayed in the plot. Currently direction of the
    principle component
    of the displayed subcomposition is displayed as a line. Later on a
    the principle componenent of the whole dataset should be displayed.}
  \item{col.pca}{The color to draw the principle component.}
  \item{logscale}{logical indication, whether logscale should be used}
  \item{xlim}{2xncol(x)-matrix giving the xlims for the columns of x}
  \item{ylim}{2xncol(x)-matrix giving the ylims for the columns of x}
}
\details{
TO DO: fix pca bug
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{plot.aplus}},
  \code{\link{qqnorm.acomp}},\code{\link{boxplot.acomp}} }

\examples{
data(SimulatedAmounts)
plot(aplus(sa.lognormals))
plot(rplus(sa.lognormals))
plot(aplus(sa.lognormals5))
plot(rplus(sa.lognormals5))
}
\keyword{multivariate}
