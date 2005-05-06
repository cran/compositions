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
          mean.acomp(x,\dots, na.action=get(getOption("na.action")))
          mean.rcomp(x,\dots, na.action=get(getOption("na.action")))
          mean.aplus(x,\dots, na.action=get(getOption("na.action")))
          mean.rplus(x,\dots, na.action=get(getOption("na.action")))
          mean.rmult(x,\dots, na.action=get(getOption("na.action")))
          }
\arguments{
  \item{x}{a classed dataset of amounts or compositions}
  \item{\dots}{further arguments to \code{\link{mean}} e.g. \code{trim}}
  \item{na.action}{ The na.action to be used: one of
    \code{\link{na.omit}},\code{\link{na.fail}},\code{\link{na.pass}}}

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

\seealso{
  \code{\link{clo}}, \code{\link{mean.col}},
  \code{\link{geometricmean}}, \code{\link{acomp}},
\code{\link{rcomp}}, \code{\link{aplus}}, \code{\link{rplus}}
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
