\name{rplus}
\alias{rplus}
\title{Amounts i.e. positive numbers analysed as objects of the real vector space}
\description{
  A class to analyse positive amounts in a classical (non-logarithmic) framework.
}
\usage{
          rplus(X,parts=1:NCOL(oneOrDataset(X)),total=NA)
          }
\arguments{
  \item{X}{vector or dataset of positive numbers considered as amounts}
  \item{parts}{vector containing the indices  xor names of the columns to be used}
  \item{total}{a numeric vectors giving the total amount of each dataset. }
}
\value{
  a vector of class \code{"rplus"} representing a vector of amounts
  or a matrix of class \code{"rplus"} representing
  multiple vectors of amounts, by rows.  
}
\details{
  Many multivariate datasets essentially describe amounts of D different
  parts in a whole. When the whole is large in relation to the
  considered parts, such that they do not exclude each other, and when
  the total amount of each componenten is actually determined by the
  phenomenon under investigation and not by sampling artifacts (such as dilution
  or sample preparation) then the parts can be treated as amounts rather
  than as a composition (cf. \code{\link{rcomp}}, \code{\link{aplus}}).
  \cr
  In principle, amounts are just real-scaled numbers with the single
  restriction that they are nonnegative. Thus they can be analysed by
  any multivariate analysis method. This class provides a simple access
  interface to do so. It tries to keep in mind the positivity
  property of amounts and the special point zero. However there are
  strong arguments why an analyis based on log-scale might be much more
  adapted to the problem. This log-approach is provided by the class
  \code{\link{aplus}}. 

  The classes rcomp, acomp, aplus, and rplus are designed in a fashion as similar as
  possible in order to allow direct comparison between results obtained  
  by the different approaches. In particular, the aplus logistic transform
  \code{\link{ilt}} is mirrored
  by the simple identity transform \code{\link{iit}}. In terms
  of computer science, this identity mapping is actually mapping an object
  of type "rplus" to a class-less datamatrix. 

  

}
\seealso{
  \code{\link{iit}},\code{\link{rcomp}}, \code{\link{aplus}},
  \code{\link{princomp.rplus}}, 
  \code{\link{plot.rplus}}, \code{\link{boxplot.rplus}},
  \code{\link{barplot.rplus}}, \code{\link{mean.rplus}},
  \code{\link{var.rplus}}, \code{\link{variation.rplus}},
  \code{\link{cov.rplus}}, \code{\link{msd}}
}

\references{

}
\examples{
data(SimulatedAmounts)
plot(rplus(sa.lognormals))

}
\keyword{multivariate}
