\name{mean.acomp}
\alias{mean.acomp}
\alias{mean.rcomp}
\alias{mean.aplus}
\alias{mean.rplus}
\alias{mean.rmult}
\title{Mean amounts and mean compositions}
\description{
  Compute the mean in the several approaches of compositional and amount
  data analysis.
}
\usage{
          \method{mean}{acomp}(x,\dots)
          \method{mean}{rcomp}(x,\dots)
          \method{mean}{aplus}(x,\dots)
          \method{mean}{rplus}(x,\dots)
          \method{mean}{rmult}(x,\dots,na.action=get(getOption("na.action")))
          }
\arguments{
  \item{x}{a classed dataset of amounts or compositions}
  \item{\dots}{further arguments to \code{\link{mean}} e.g. \code{trim}}
  \item{na.action}{na.action }
}
\value{
  The mean is given as a composition or amount vector of the same class as the original dataset.
}
\details{
  The different compositional approaches \code{\link{acomp}},
  \code{\link{rcomp}},
  \code{\link{aplus}}, \code{\link{rplus}} correpond to different
  geometries. The mean is calculated in the respective canonical
  geometry by applying a canonical transform (see \code{\link{cdt}}), taking ordinary
  \code{\link{mean.col}} and backtransforming.\cr

  The Aitchison geometries imply that \code{mean.acomp} and \code{mean.aplus} are
  geometric means, the first one closed. The real geometry implies that
 \code{mean.rcomp} and \code{mean.rplus} are arithmetic means, the first
  one resulting in a closed composition.\cr

In all cases the mean is again an object of the same class. 


}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{clo}}, \code{\link{mean.col}},
  \code{\link{geometricmean}}, \code{\link{acomp}},
\code{\link{rcomp}}, \code{\link{aplus}}, \code{\link{rplus}}
}
\section{Missing Policy}{
  For the additive scales (rcomp,rplus) the SZ and BDT are 
  treated as zeros and MAR and MNAR as missing information. 
  This is not strictly correct for MNAR. 
  \cr
  For relative scales (acomp,aplus), all four types of missings 
  are treated as missing information. This corresponds to the 
  idea that BDT are truncated values (and have the correspoding 
  effect in taking means). For SZ and MAR, only the components in 
  the observed subcomposition are fully relevant. Finally, for MNAR 
  the problem is again that nothing could be done without knowing
  the MNAR mechanism, so the analysis is limited to taking them as 
  MAR, and being \emph{careful} with the interpretation.
  Missing and Below Detecion Limit Policy is explained in more detail 
  in \link{compositions.missing}. 
  }
\examples{
data(SimulatedAmounts)
mean.col(sa.lognormals)
mean(acomp(sa.lognormals))
mean(rcomp(sa.lognormals))
mean(aplus(sa.lognormals))
mean(rplus(sa.lognormals))
mean(rmult(sa.lognormals))
}
\keyword{multivariate}
