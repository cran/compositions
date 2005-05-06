\name{meanrow}
\alias{mean.row}
\alias{mean.col}

\title{The arithmetic mean of rows or columns}
\description{
  Computes the arithmetic mean.
}
\usage{
          mean.row(x,\dots, na.action=get(getOption("na.action")))
          mean.col(x,\dots, na.action=get(getOption("na.action")))
          }
\arguments{
  \item{x}{a numeric vector or matrix of data }
  \item{\dots}{arguments to \code{\link{mean}} }
  \item{na.action}{ The na.action to be used: one of
    \code{\link{na.omit}},\code{\link{na.fail}},\code{\link{na.pass}}}
}
\value{
  The arithmetic means of the rows
  (mean.row)
  or  columns (mean.col) of x.
}
\details{

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{mean.rplus}} }

\examples{
data(SimulatedAmounts)
mean.col(sa.tnormals)
mean.row(sa.tnormals)

}
\keyword{univar}
