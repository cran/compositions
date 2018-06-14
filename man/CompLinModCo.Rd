\name{CompLinModCoReg}
\alias{CompLinModCoReg}
\title{Compositional Linear Model of Coregionalisation}
\description{
  Creates a Variogram model according to the linear model of spatial
  corregionalisation for a compositional geostatistical analysis.
}
\usage{
CompLinModCoReg(formula,comp,D=ncol(comp),envir=environment(formula))
          }
\arguments{
  \item{formula}{A formula without left side providing a formal
    description of a variogram model.}
  \item{comp}{a compositional dataset, needed to provide the frame size}
  \item{D}{The dimension of the multivariate dataset}
  \item{envir}{The enviroment in which formula should be interpreted.}
}
\value{
  A variogram function. 
}
\details{
  The linear model of coregionalisation says uses the fact, that sums of
  valid variogram models are valid variograms and that skalar variograms
  multiplid with a positiv definit matrix are valid vector
  variograms. The function computes such a variogram function from a
  formal description. The formula is a sum. Each summand is either a
  product of a matrix description and a scalar variogram description or
  only a scalar variogram description. Scalar variogram descriptions are
  either formal function calls to
  \itemize{
    \item{\code{sph(range)}}{ for sperical variogram with range \code{range}}
    \item{\code{exp(range)}}{ for an exponential variogram with range
      \code{range}}
    \item{\code{gauss(range)}}{ for a Gaussian variogram with range
    \code{range}}
    \item{\code{gauss(range)}}{ for a cardinal sine variogram with (non-effective) range
    \code{range}}
    \item{\code{pow(range)}}{ for an power variogram with range
      \code{range}}
    \item{\code{lin(unit)}}{ linear variogram 1 at \code{unit}.}
    \item{\code{nugget()}}{ for adding a nuggeteffect.}
  }
  alternatively it can be any expression, which will be evaluated in
  envir and should depende on a dataset of distantce vectrs \code{h}.
  \cr
  The matrix description always comes first. It can be \code{R1} for a
  rank 1 matrix \code{PSD} for a Positive Semidefinite matrix or \(S\)
  for a scalar Sill factor to be multipled with identiy. Or any other
  construct evaluating to a matrix or a function of some parameters with
  default values, that if called is evaluated to a positive semidefinit
  matrix. \code{R1} and \code{PSD} can also be written as calls with
  providing a vector or respectively a matrix providing the parameter.
  \cr
  The variogram is created with default parametervalues. The parameters
  can later be modified by modifiying the default parameter with
  assignments like \code{formals(vg)$sPSD1 =
    parameterPosdefMat(4*diag(5))}.
  We would anyway expect you to fit the model to the data by a command
  like \code{vgmFit(logratioVariogram(...),CompLinModCoReg(...))}  
}
\references{
  What to cite??
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{vgram2lrvgram}},
  \code{\link{CompLinModCoReg}},
  \code{\link{vgmFit}}
}
\examples{
\dontrun{
data(juraset)
X <- with(juraset,cbind(X,Y))
comp <- acomp(juraset,c("Cd","Cu","Pb","Co","Cr"))
CompLinModCoReg(~nugget()+sph(0.5)+R1*exp(0.7),comp)
CompLinModCoReg(~nugget()+R1*sph(0.5)+R1*exp(0.7)+(0.3*diag(5))*gauss(0.3),comp)
CompLinModCoReg(~nugget()+R1*sph(0.5)+R1(c(1,2,3,4,5))*exp(0.7),comp)
}
}
\keyword{multivariate}
