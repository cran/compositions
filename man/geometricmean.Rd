\name{geometricmean}
\alias{geometricmean}
\alias{geometricmean.row}
\alias{geometricmean.col}
\alias{gsi.geometricmean}
\alias{gsi.geometricmean.row}
\alias{gsi.geometricmean.col}
\title{The geometric mean}
\description{
  Computes the geometric mean.
}
\usage{
          geometricmean(x,\dots)
          geometricmean.row(x,\dots)
          geometricmean.col(x,\dots)
          gsi.geometricmean(x,\dots)
          gsi.geometricmean.row(x,\dots)
          gsi.geometricmean.col(x,\dots)
          }
\arguments{
  \item{x}{a numeric vector or matrix of data }
  \item{\dots}{further arguments to compute the mean }
}
\value{
  The geometric means of x as a whole (geometricmean), its rows
  (geometricmean.row) or its columns (geometricmean.col).
}
\details{
  The geometric mean is defined as:
  \deqn{geometricmean(x) := \left( \prod_{i=1}^n x_i\right)^{1/n}}{
    geometricmean(x) := \code{prod(x)^(1/length(x))}
  }
  The geometric mean is actually computed by
  \code{exp(mean(log(c(unclass(x))),...))}.
}
\section{Missing Policy}{
  The the first three functions take the geometric mean of all non-missing values. 
  This is because they should yield a result in term of data analysis. 
 
  Contrarily, the gsi.* functions inherit the arithmetic IEEE policy of R through
  \code{exp(mean(log(c(unclass(x))),...))}. Thus,  NA codes a not available i.e. 
  not measured, NaN codes a below detection limit, and 0.0 codes a structural zero.
  If any of the elements involved is 0, NA or NaN the result is of the same
  type. Here 0 takes precedence over NA, and NA takes precedence
  over NaN. For example, if a structural 0 appears, the geometric mean is 0
  regardless of the presence of NaN's or NA's in the rest. Values below detection 
  limit become NaN's if they are coded as negative values.
}

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{mean.rplus}}
}
\examples{
geometricmean(1:10)  
geometricmean(c(1,0,NA,NaN))  # 0
X <- matrix(c(1,NA,NaN,0,1,2,3,4),nrow=4)
X  
geometricmean.row(X)
geometricmean.col(X)
}
\keyword{univar}
