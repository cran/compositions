\name{clr2ilr}
\alias{clr2ilr}
\alias{ilr2clr}
\alias{clrvar2ilr}
\alias{ilrvar2clr}
\title{Convert between clr and ilr, and between cpt and ipt. Acts in vectors and in bilinear forms.}
\description{
  Compute the centered log ratio transform of a (dataset of)
  isometric log-ratio transform(s) and its inverse. Equivalently, 
  compute centered and isometric planar transforms from each other.
}
\usage{
clr2ilr( x , V=ilrBase(x) )
ilr2clr( z , V=ilrBase(z=z) )
clrvar2ilr( varx , V=ilrBase(D=ncol(varx)) )
ilrvar2clr( varz , V=ilrBase(D=ncol(varz)+1) )

          }
\arguments{
  \item{x}{the clr/cpt-transform of composition(s)}
  \item{z}{the ilr/ipt-transform of composition(s)}
  \item{varx}{variance or covariance matrix of clr/cpt-transformed
    compositions}
  \item{varz}{variance or covariance matrix of ilr/ipt-transformed
    compositions}
  \item{V}{a matrix with columns giving the chosen  basis of the clr-plane}
}
\value{
  \code{clr2ilr} gives the ilr/ipt transform of the same composition(s),\cr
  \code{ilr2clr} gives the clr/cpt transform of the same
  composition(s),\cr
  \code{clrvar2ilr} gives the variance-/covariance-matrix of the ilr/ipt transform of the same compositional data set,\cr
  \code{ilrvar2clr} gives the variance-/covariance-matrix of the clr/cpt
  transform of the same compositional data set.\cr
  
}
\details{
  These functions perform a matrix multiplication with \code{V} in an
  appropriate way. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{ilr}}, \code{\link{ipt}}
}
\references{
  Egozcue J.J., V. Pawlowsky-Glahn, G. Mateu-Figueras and
  C. Barcel'o-Vidal (2003) Isometric logratio transformations for
  compositional data analysis. \emph{Mathematical Geology}, \bold{35}(3)
  279-300\cr
  Aitchison, J, C. Barcel'o-Vidal, J.J. Egozcue, V. Pawlowsky-Glahn
  (2002) A consise guide to the algebraic geometric structure of the
  simplex, the sample space for compositional data analysis, \emph{Terra
    Nostra}, Schriften der Alfred Wegener-Stiftung, 03/2003
  

}
\examples{
data(SimulatedAmounts)
ilr.inv(clr2ilr(clr(sa.lognormals)))-clo(sa.lognormals)
clr.inv(ilr2clr(ilr(sa.lognormals)))-clo(sa.lognormals)
ilrvar2clr(var(ilr(sa.lognormals)))-var(clr(sa.lognormals))
clrvar2ilr(var(cpt(sa.lognormals)))-var(ipt(sa.lognormals))

}
\keyword{multivariate}
