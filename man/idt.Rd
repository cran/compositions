\name{idt}
\alias{idt}
\alias{idt.default}
\alias{idt.acomp}
\alias{idt.rcomp}
\alias{idt.aplus}
\alias{idt.rplus}
\alias{idt.rmult}
\alias{idt.factor}
\alias{idt.inv}
\alias{idt.inv.default}
\alias{idt.inv.acomp}
\alias{idt.inv.rcomp}
\alias{idt.inv.aplus}
\alias{idt.inv.rplus}
\alias{idt.inv.rmult}
\title{Isometric default transform}
\description{
  Compute the isometric default transform of a vector (or dataset) of
  compositions or amounts in the selected class.
}
\usage{
          idt(x,...)
          \method{idt}{default}( x,... )
          \method{idt}{acomp}( x ,...)
          \method{idt}{rcomp}( x ,...)
          \method{idt}{aplus}( x ,...)
          \method{idt}{rplus}( x ,...)
          \method{idt}{rmult}( x ,...)
          \method{idt}{factor}( x ,...)
          idt.inv(x,orig,...)
          \method{idt.inv}{default}( x ,orig,...)
          \method{idt.inv}{acomp}( x ,orig,...)
          \method{idt.inv}{rcomp}( x ,orig,...)
          \method{idt.inv}{aplus}( x ,orig,...)
          \method{idt.inv}{rplus}( x ,orig,...)
          \method{idt.inv}{rmult}( x ,orig,...)

          }
\arguments{
  \item{x}{a classed amount or composition, to be transformed with its
    isometric default transform, or its inverse}
  \item{\dots}{generic arguments past to underlying functions}
  \item{orig}{a compositional object which should be mimicked
             by the inverse transformation. It is used to determine the
	     backtransform to be used, and eventually to
	     reconstruct the names of the parts. It is the generic
	     argument. Typically the orig-argument is the dataset that
	     has been transformed in the first place.}

}
\value{
  A corresponding matrix of row-vectors containing the transforms.
}
\details{
  The general idea of this package is to analyse the same data with
  different geometric concepts, in a fashion as similar as possible. For each of the
  four concepts there exists an isometric transform expressing the geometry
  in a full-rank euclidean vector space. Such a transformation is computed 
  by \code{idt}. For \code{\link{acomp}} the transform is \code{\link{ilr}}, for
  \code{\link{rcomp}} it is \code{\link{ipt}}, for 
  \code{\link{aplus}} it is \code{\link{ilt}}, and for
  \code{\link{rplus}} it is \code{\link{iit}}. Keep in mind that the
  transform does not keep the variable names, since there is no guaranteed 
  one-to-one relation between the original parts and each transformed
  variable.
  \cr
  The inverse \code{idt.inv} is intended to allow for an "easy" and automatic 
  back-transformation, without intervention of the user. The argument \code{orig} 
  (the one determining the behaviour of \code{idt.inv} as a generic function) 
  tells the function which back-transformation should be applied, and
  gives the column names of \code{orig} to the back-transformed
  values of \code{x}. Therefore, it is very conventient to give the original classed
  data set used in the analysis as \code{orig}.
}
\author{R. Tolosana-Delgado, K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{cdt}}, \code{\link{ilr}}, \code{\link{ipt}},
  \code{\link{ilt}}, \code{\link{cdt.inv}}, \code{\link{ilr.inv}}, \code{\link{ipt.inv}}, 
            \code{\link{ilt.inv}}, \code{\link{iit.inv}}
}
\references{
    van den Boogaart, K.G. and R. Tolosana-Delgado (2007) "compositions": a unified 
  R package to analyze Compositional Data, \emph{Computers & Geosciences}. 
  (in press).
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
  data(Hydrochem)
  x = Hydrochem[,c("Na","K","Mg","Ca")]
  y = acomp(x)
  z = idt(y)
  y2 = idt.inv(z,y)
  par(mfrow=c(2,2))
  for(i in 1:4){plot(y[,i],y2[,i])}
}
\keyword{multivariate}
