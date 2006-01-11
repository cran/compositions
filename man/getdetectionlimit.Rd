\name{getdetectionlimit}
\alias{getDetectionlimit}
\title{Gets the detection limit stored in the data set}
\description{
  The detection limit of those values below-detection-limit are 
  stored as negative values in compositional dataset. This function 
  extracts that information.
}
\usage{
getDetectionlimit(x,dl=attr(x,"detectionlimit"))
}
\arguments{
  \item{x}{a data set}
  \item{dl}{a default to replace the information in the dataset}
}
\value{
  a matrix in the same shape as x, with a positive value (the detection limit)
  where available, and NA in the other cells. 
}
\details{
  For a proper treatment of truncated data it would be necessary to know
  the detection limit even for observed data. Unfortunately, there
  is no clear way to encode this information without annoying the user.
}
\references{
  
}
\seealso{
  \code{\link{compositions.missings}},\code{\link{zeroreplace}}
}
\examples{
x <- c(2,-0.5,4,3,-0.5,5,BDLvalue,MARvalue,MNARvalue)
getDetectionlimit(x)

}
\author{}
\keyword{}