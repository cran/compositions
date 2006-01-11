\name{simulatemissings}
\alias{simulateMissings}
\title{Artifical simulation of various kinds of missings}
\description{
  This is a simulation mechanism to check the missing technique. It just
  generates additional missings of the various types in a given dataset.
}
\usage{
simulateMissings(x,detectionlimit=NULL,
     knownlimit=FALSE,MARprob=0.0,MNARprob=0.0,
     mnarity=0.5,SZprob=0.0)
}
\arguments{
  \item{x}{a dataset that should get the missings}
  \item{detectionlimit}{the detection limit described in
    \code{\link{clo}}, to impose an artificial detection limit}
  \item{knownlimit}{a boolean indicating wether the actual detection
    limit is still known in the dataset.}
  \item{MARprob}{the probability of occurence of 'Missings At Random' values}
  \item{MNARprob}{the probability of occurrence of 'Missings Not At Random'. 
    The tendency is that small values have a higher probability to
    be missed.}
  \item{mnarity}{a number between 0 and 1 giving the strength of the
    influence of the actual value in becoming a MNAR. 0 means a MAR
    like behavior and 1 means that it is just the smallest values that
    is lost}
  \item{SZprob}{the probability to obtain a structural zero. This is done
    at random like a MAR.}
}
\value{
  A dataset like \code{x} but with some additional missings.

}
\details{
  Without any additional parameters no missings are generated. 
  The procedure to generate MNAR affects all variables.
}
\references{
  See \link{compositions.missings} for more details.
}
\seealso{
  \link{compositions.missings}
}
\examples{
data(SimulatedAmounts)
x <- acomp(sa.lognormals)
xnew <- simulateMissings(x,detectionlimit=0.05,MAR=0.05,MNAR=0.05,SZ=0.05)
xnew
plot(missingSummary(xnew))
}
\author{K.Gerald van den Boogaart}
\keyword{datagen}