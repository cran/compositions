\name{cdt}
\alias{cdt}
\alias{cdt.default}
\alias{cdt.acomp}
\alias{cdt.rcomp}
\alias{cdt.aplus}
\alias{cdt.rplus}
\alias{cdt.rmult}
\alias{cdt.factor}
\title{Centered default transform}
\description{
  Compute the centered default transform of a (dataset of)
  composition or amount.
}
\usage{
          cdt(x)
          \method{cdt}{default}( x )
          \method{cdt}{acomp}( x )
          \method{cdt}{rcomp}( x )
          \method{cdt}{aplus}( x )
          \method{cdt}{rplus}( x )
          \method{cdt}{rmult}( x )
          \method{cdt}{factor}( x )
          }
\arguments{
  \item{x}{a classed amount or composition (or a matrix of), to be transformed with its
    centered default transform}
}
\value{
  A corresponding matrix or vector containing the transforms.
}
\details{
  The general idea of this package is to analyse the same data with
  different geometric concepts as similar as possible. For each of the
  four concepts there exists a unique transform expressing the geometry
  in a linear subspace, keeping the relation to the variables. This
  unique transformation is computed by \code{cdt}. For
  \code{\link{acomp}} the transform is \code{\link{clr}}, for
\code{\link{rcomp}} it is \code{\link{cpt}}, for 
  \code{\link{aplus}} it is \code{\link{ilt}}, and for
  \code{\link{rplus}} it is \code{\link{iit}}. Each component of the result
  is identified with a unit vector in the direction of the corresponding
  component of the original composition or amount. Keep in mind that the
  transform is not necessarily surjective and thus variances in the
  image space might be singular.

  
}
\seealso{
  \code{\link{idt}}, \code{\link{clr}}, \code{\link{cpt}}, \code{\link{ilt}}, \code{\link{iit}}
}
\references{
  

}
\examples{
\dontrun{
# the cdt is defined by
cdt         <- function(x) UseMethod("cdt",x)
cdt.default <- function(x) x
cdt.acomp   <- clr 
cdt.rcomp   <- cpt 
cdt.aplus   <- ilt 
cdt.rplus   <- iit 
}
cdt(acomp(1:5))
cdt(rcomp(1:5))

}
\keyword{multivariate}
