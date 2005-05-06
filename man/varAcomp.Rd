\name{var.acomp}
\alias{var}
\alias{var.default}
\alias{var.acomp}
\alias{var.rcomp}
\alias{var.aplus}
\alias{var.rplus}
\alias{var.rmult}
\alias{cov}
\alias{cov.default}
\alias{cov.acomp}
\alias{cov.rcomp}
\alias{cov.aplus}
\alias{cov.rplus}
\alias{cov.rmult}
\title{Variances and covariances of amounts and compositions}
\description{
  Compute the (co)variance matrix in the several approaches of compositional
 and amount data analysis.
}
\usage{
          var(x,\dots)
          \method{var}{default}(x, y=NULL, na.rm=FALSE, use, \dots)
          \method{var}{acomp}(x,y=NULL,\dots)
          \method{var}{rcomp}(x,y=NULL,\dots)
          \method{var}{aplus}(x,y=NULL,\dots)
          \method{var}{rplus}(x,y=NULL,\dots)
          \method{var}{rmult}(x,y=NULL,\dots)
          cov(x,y=x,\dots)
          \method{cov}{default}(x, y=NULL, use="all.obs", method=c("pearson",
    "kendall", "spearman"), \dots)
          \method{cov}{acomp}(x,y=NULL,\dots)
          \method{cov}{rcomp}(x,y=NULL,\dots)
          \method{cov}{aplus}(x,y=NULL,\dots)
          \method{cov}{rplus}(x,y=NULL,\dots)
          \method{cov}{rmult}(x,y=NULL,\dots)
          }
\arguments{
  \item{x}{a dataset, eventually of amounts or compositions}
  \item{y}{a second dataset, eventually of amounts or compositions}
  \item{na.rm}{see \code{\link[stats]{var}}}
  \item{use}{see \code{\link[stats]{var}}}
  \item{method}{see \code{\link[stats]{cov}}}
  \item{\dots}{further arguments to \code{\link[stats]{var}}
    e.g. \code{use}}
}
\value{
  The variance matrix of x or the covariance matrix of x and
  y.}
\details{
  The basic functions of \code{\link[stats]{var}},
  \code{\link[stats]{cov}} are turned to
  S3-generics. The original versions are copied to the default
  method. This allows us to introduce generic methods to handle
  variances and covariances of other datatypes such as amounts or
  compositions.
  \cr

  If classed amounts or compositions are involved, they are transformed
  with their corresponding transforms, using the centered default
  transform (\code{\link{cdt}}). That implies that the variances have to
  be interpreded in a log scale level for  \code{\link{acomp}} and
  \code{\link{aplus}}.
  \cr
  We should be aware that variance matrices of compositions are
  singular. They can be transformed to the correponding nonsingular
  variances of ilr or ipt -space by \code{\link{clrvar2ilr}}.

  In R versions older than v2.0.0, \code{\link[base]{var}} and
  \code{\link[base]{cov}} were defined in package ``base'' instead of in ``stats''.
  This might produce some misfunction.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{cdt}}, \code{\link{clrvar2ilr}}, \code{\link{clo}},
\code{\link{mean.acomp}}, \code{\link{acomp}}, \code{\link{rcomp}},
\code{\link{aplus}}, \code{\link{rplus}}, \code{\link{variation}}
}

\examples{
data(SimulatedAmounts)
mean.col(sa.lognormals)
var(acomp(sa.lognormals))
var(rcomp(sa.lognormals))
var(aplus(sa.lognormals))
var(rplus(sa.lognormals))
cov(acomp(sa.lognormals5[,1:3]),acomp(sa.lognormals5[,4:5]))
cov(rcomp(sa.lognormals5[,1:3]),rcomp(sa.lognormals5[,4:5]))
cov(aplus(sa.lognormals5[,1:3]),aplus(sa.lognormals5[,4:5]))
cov(rplus(sa.lognormals5[,1:3]),rplus(sa.lognormals5[,4:5]))
cov(acomp(sa.lognormals5[,1:3]),aplus(sa.lognormals5[,4:5]))

svd(var(acomp(sa.lognormals)))

}
\keyword{multivariate}
