\name{cor.acomp}
\alias{cor}
\alias{cor.default}
\alias{cor.acomp}
\alias{cor.rcomp}
\alias{cor.aplus}
\alias{cor.rplus}
\alias{cor.rmult}
\title{Correlations of amounts and compositions}
\description{
  Computes the correlation matrix in the various approaches of compositional
 and amount  data analysis.
}
\usage{
          cor(x,y=NULL,\dots)
          \method{cor}{default}(x, y=NULL, use="all.obs", method=c("pearson", 
    "kendall", "spearman"),\dots)
          \method{cor}{acomp}(x,y=NULL,\dots)
          \method{cor}{rcomp}(x,y=NULL,\dots)
          \method{cor}{aplus}(x,y=NULL,\dots)
          \method{cor}{rplus}(x,y=NULL,\dots)
          \method{cor}{rmult}(x,y=NULL,\dots)
          }
\arguments{
  \item{x}{a dataset, eventually of amounts or compositions}
  \item{y}{a second dataset, eventually of amounts or compositions}
  \item{use}{see \code{\link[stats]{cor}}}
  \item{method}{see \code{\link[stats]{cor}}}
  \item{\dots}{further arguments to \code{\link[stats]{cor}}
    e.g. \code{use}}
}
\value{
  The correlation  matrix.
}
\details{
The correlation matrix does not make much sense for compositions.

  In R versions older than v2.0.0, \code{\link[base]{cor}} was defined
  in package ``base'' instead of in ``stats''.
  This might produce some misfunction.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{var.acomp}}
}

\examples{
data(SimulatedAmounts)
mean.col(sa.lognormals)
cor(acomp(sa.lognormals5[,1:3]),acomp(sa.lognormals5[,4:5]))
cor(rcomp(sa.lognormals5[,1:3]),rcomp(sa.lognormals5[,4:5]))
cor(aplus(sa.lognormals5[,1:3]),aplus(sa.lognormals5[,4:5]))
cor(rplus(sa.lognormals5[,1:3]),rplus(sa.lognormals5[,4:5]))
cor(acomp(sa.lognormals5[,1:3]),aplus(sa.lognormals5[,4:5]))


}
\keyword{multivariate}
