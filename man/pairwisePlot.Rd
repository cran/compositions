\name{pairwiseplot}
\alias{pairwisePlot}
\alias{pairwisePlot.default}
\title{Creates a paneled plot like pairs for two different datasets. }
\description{
  Creates a plot for each element of two lists or each column of each
  dataset against each of the second.
}
\usage{
pairwisePlot(X,Y,...)
\method{pairwisePlot}{default}(X,Y=X,...,xlab=deparse(substitute(X)),ylab=deparse(substitute(Y)),nm=c(length(Y),length(X)),panel=plot)
}
\arguments{
  \item{X}{a list, a data.frame, or a matrix representing the first set
          of things to be displayed.}
  \item{Y}{a list, a data.frame, or a matrix representing the second set
          of things to be displayed.}
  \item{\dots}{furter parameters to the panel function}
  \item{xlab}{The sequence of labels for the elements of
          X. Alternatively the labels can be given as colnames or names
          of X. This option takes precedence if specified.}
  \item{ylab}{The sequence of labels for the elements of
          Y. Alternatively the labels can be given as colnames or names
          of Y. This option takes precedence if specified.}
  \item{nm}{the parameter to be used in the call
           \code{par(mfrow=nm)}. If NULL no parameter is setted and a
       sequence of plots can be generated. }
  \item{panel}{The panel function to plot the individual panels.
    it is called as \code{panel(x,y,xlab=xlab,ylab=ylab,...)}. Thus the
    panel function must be capable of of taking these arguments. It must
     also set up its own plot. There is no negotiation on coordinate
     system. 
   }
  }
\details{ This is a light weight conveniance function to several
aspects of one dataset against several aspects of another dataset.  It
is far more straight forward than e.g. the \code{pairs} function and
does not do any internal computation rather than organizing the names.
}
\references{
  Boogaart, K.G. v.d. , R. Tolosana (2008) Mixing Compositions and Other
  scales, Proceedings of CodaWork 08, to appear

  
  \url{http://ima.udg.es/Activitats/CoDaWork03}\cr
  
  \url{http://ima.udg.es/Activitats/CoDaWork05}

  \url{http://ima.udg.es/Activitats/CoDaWork08}

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}

\seealso{
  \code{\link{plot.aplus}}, \code{\link{balance}}
}
\examples{

X <- rnorm(100)
Y <- rnorm.acomp(100,acomp(c(A=1,B=1,C=1)),0.1*diag(3))+acomp(t(outer(c(0.2,0.3,0.4),X,"^")))

pairs(cbind(ilr(Y),X),panel=function(x,y,...) {points(x,y,...);abline(lm(y~x))})
pairs(cbind(balance(Y,~A/B/C),X),panel=function(x,y,...) {points(x,y,...);abline(lm(y~x))})
pairwisePlot(balance(Y,~A/B/C),X)
pairwisePlot(X,balance(Y,~A/B/C),panel=function(x,y,...) {plot(x,y,...);abline(lm(y~x))})
pairwisePlot(X,balance01(Y,~A/B/C))

# A function to extract a portion representation of subcompsitions
# with two elements:
subComps <- function(X,...,all=list(...)) {
  X <- oneOrDataset(X)
  nams <- sapply(all,function(x) paste(x[[2]],x[[3]],sep=","))
  val  <- sapply(all,function(x) X[,match(as.character(x[[2]]),colnames(X)) ]/(X[,match(as.character(x[[2]]),colnames(X)) ]+X[,match(as.character(x[[3]]),colnames(X)) ]))
  colnames(val)<-nams
  val
}

pairwisePlot(X,subComps(Y,A~B,A~C,B~C))

}
\keyword{hplot}
