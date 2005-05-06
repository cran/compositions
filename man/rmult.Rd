\name{rmult}
\alias{rmult}
\title{Simple treatment of real vectors }
\description{
  A class to analyse real multivariate vectors.
}
\usage{
          rmult(X,parts=1:NCOL(oneOrDataset(X)))
          }
\arguments{
  \item{X}{vector or dataset of numbers considered as elements of a R-vector}
  \item{parts}{vector containing the indices  xor names of the columns to be used}
}
\value{
  a vector of class \code{"rmult"} representing one vector
  or a matrix of class \code{"rmult"}, representing
  multiple vectors by rows.  
}
\details{
The \code{rmult} class is a simple convenience class to treat
data in the scale of real vectors just like data in the scale of real
numbers. A major aspect to take into account is that the internal arithmetic of R is
different for these vectors.    
}
\seealso{
  \code{\link{+.rmult}}, \code{\link{scalar}}, \code{\link{norm}},
  \code{\link{\%*\%.rmult}},
  \code{\link{rplus}}, \code{\link{acomp}},
}

\examples{
plot(rnorm.rmult(30,mean=0:4,var=diag(1:5)+10))

}
\keyword{multivariate}
