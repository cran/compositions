\name{endpointCoordinates}
\alias{endpointCoordinates}
\alias{endpointCoordinates.inv}
\alias{endpointCoordinates.default}
\alias{endpointCoordinates.acomp}
\alias{endpointCoordinates.aplus}
\alias{endpointCoordinates.rplus}
\alias{endpointCoordinatesInv}
\alias{endpointCoordinatesInv.rmult}
\alias{endpointCoordinatesInv.acomp}
\alias{endpointCoordinatesInv.aplus}
\alias{endpointCoordinatesInv.rcomp}
\alias{endpointCoordinatesInv.rplus}
\title{Amounts in barytic-coordinates}
\description{
  Computes the convex combination of amounts given by \code{endpoints}
  to explain \code{X} as good as possible.
}
\usage{

          endpointCoordinates(X,\dots)
          endpointCoordinatesInv(K,endpoints,...)
          \method{endpointCoordinates}{default}(X,endpoints=diag(gsi.getD(X)), \dots)
          \method{endpointCoordinates}{acomp}(X,endpoints=clr.inv(diag(gsi.getD(X))),...)
          \method{endpointCoordinates}{aplus}(X,endpoints,...)
          \method{endpointCoordinates}{rplus}(X,endpoints,...)
          \method{endpointCoordinatesInv}{rmult}(K,endpoints,...)
          \method{endpointCoordinatesInv}{acomp}(K,endpoints,...)
          \method{endpointCoordinatesInv}{rcomp}(K,endpoints,...)
          \method{endpointCoordinatesInv}{aplus}(K,endpoints,...)
          \method{endpointCoordinatesInv}{rplus}(K,endpoints,...)
          }
\arguments{
  \item{X}{a dataset of amounts or compositions, to be represented in
    as convex combination of the endpoints in the given geometry}
  \item{K}{Konvex combination weights to the endpoints}
  \item{endpoints}{a dataset of extremal compositions from the same
    space as X. The number of endpoints given must not exceed the
    dimension of the space plus one.}
  \item{\dots}{currently unused}
}
\value{
  The \code{endpointCoordinates} functions give a \code{"rmult"}-dataset
  giving the convex weights, which allow to combine \code{X} from
  \code{endpoints} as good as possible. The result is an \code{"rmult"}
  since there is guarantee that the resulting weights are positive.
  \cr
  The \code{endpointCoordinates} functions reconstruct the convex
  combination from coordinates \code{K} and the given
  \code{endpoints}. The class of \code{endpoints} determines the
  geometry chosen and the class of the result.  
}
\details{
  The convex combination is performed in the respective geometry. This
  means that for rcomp positivity of the result is only guaranteed with
  extermal endmembers and that in acomp-geometry it is not possible to
  give extremal endmembers.  
  \cr
  The main idea behind this functions is that the actually observed
  composition came from a convex combination of some extremal
  compositions specified by endpoints. Strictly speaking this is
  meaningfull in strictly this sense only in rplus-geometry and under
  some special circumstances in rcomp geometry. It is not
  meaningfull in terms of mass conservation in acomp- and aplus-geometry
  due to the non mass-balancing
  character of the geometry. In rcomp-geometry it dependent on unit of
  measurements and different for volume and mass \% and only valid if
  the whole composition is observed. 
  

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de},
  Raimon Tolosana-Delgado}

\references{
  Shurtz, Robert F., 2003. Compositional geometry and mass conservation.
  Mathematical Geology 35~(8), 972--937.
}
\examples{
data(SimulatedAmounts)
ep <- aplus(rbind(c(2,1,2),c(2,2,1),c(1,2,2)))
dat <- endpointCoordinatesInv(acomp(sa.lognormals),acomp(ep))
plot(dat)
plot( acomp(endpointCoordinates(dat,acomp(ep))))

dat <- endpointCoordinatesInv(rcomp(sa.lognormals),rcomp(ep))
plot(dat)
plot( rcomp(endpointCoordinates(dat,rcomp(ep))))

dat <- endpointCoordinatesInv(aplus(sa.lognormals),aplus(ep))
plot(dat)
plot( endpointCoordinates(dat,aplus(ep)))

dat <- endpointCoordinatesInv(rplus(sa.lognormals),rplus(ep))
plot(dat)
plot(endpointCoordinates(rplus(dat),rplus(ep)))


}
\keyword{multivariate}
