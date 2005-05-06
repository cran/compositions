\name{idt}
\alias{idt}
\alias{idt.default}
\alias{idt.acomp}
\alias{idt.rcomp}
\alias{idt.aplus}
\alias{idt.rplus}
\alias{idt.rmult}
\alias{idt.factor}
\title{Isometric default transform}
\description{
  Compute the isometric default transform of a vector (or dataset) of
  compositions or amounts in the selected class.
}
\usage{
          idt(x)
          \method{idt}{default}( x )
          \method{idt}{acomp}( x )
          \method{idt}{rcomp}( x )
          \method{idt}{aplus}( x )
          \method{idt}{rplus}( x )
          \method{idt}{rmult}( x )
          \method{idt}{factor}( x )
          }
\arguments{
  \item{x}{a classed amount or composition, to be transformed with its
    isometric default transform}
}
\value{
  A corresponding matrix of row-vectors containing the transforms.
}
\details{
  The general idea of this package is to analyse the same data with
  different geometric concepts, in a fashion as similar as possible. For each of the
  four concepts there exists an isometric transform expressing the geometry
  in a full-rank euclidean vector space. Such a transformation is computed by \code{idt}. For
  \code{\link{acomp}} the transform is \code{\link{ilr}}, for
\code{\link{rcomp}} it is \code{\link{ipt}}, for 
  \code{\link{aplus}} it is \code{\link{ilt}}, and for
  \code{\link{rplus}} it is \code{\link{iit}}. Keep in mind that the
  transform does not keep the variable names, since there is no guaranteed 
  one-to-one relation between the original parts and each transformed variable.
  
}
\seealso{
  \code{\link{cdt}}, \code{\link{ilr}}, \code{\link{ipt}}, \code{\link{ilt}}, \code{\link{iit}}
}
\references{
  

}
\examples{
\dontrun{
# the idt is defined by
idt         <- function(x) UseMethod("idt",x)
idt.default <- function(x) x
idt.acomp   <- function(x) ilr(x) 
idt.rcomp   <- function(x) ipt(x) 
idt.aplus   <- ilt 
idt.rplus   <- iit 
}
idt(acomp(1:5))
idt(rcomp(1:5))

}
\keyword{multivariate}
