\name{isoPortionLines}
\alias{isoPortionLines}
\alias{isoPortionLines.acomp}
\alias{isoPortionLines.rcomp}
\alias{isoProportionLines}
\alias{isoProportionLines.acomp}
\alias{isoProportionLines.rcomp}
\title{Isoportion- and Isoproportion-lines }
\description{
  Add lines of equal portion and proportion to ternary diagrams.
}
\usage{
isoPortionLines(...)
isoPortionLines.acomp(by=0.2,at=seq(0,1,by=by),...,parts=1:3,total=1,labs=TRUE,lines=TRUE,unit="")
isoPortionLines.rcomp(by=0.2,at=seq(0,1,by=by),...,parts=1:3,total=1,labs=TRUE,lines=TRUE,unit="")
isoProportionLines(...)
isoProportionLines.acomp(by=0.2,at=seq(0,1,by=by),...,parts=1:3,labs=TRUE,lines=TRUE) 
isoProportionLines.rcomp(by=0.2,at=seq(0,1,by=by),...,parts=1:3,labs=TRUE,lines=TRUE) 
          }
\arguments{
  \item{\dots}{graphical arguments}
  \item{at}{numeric in [0,1]: Which portions/proportions should be marked?}
  \item{by}{numeric in (0,1]:Steps between protions/proportions}
  \item{parts}{numeric vector subset of \{1,2,3\}: the variables to be
    marked}
  \item{total}{The total amount to be used in labeling}
  \item{labs}{logical: plot the labels?}
  \item{lines}{logical: plot the lines?}
  \item{unit}{marking the units e.g. "\%"}

}
\details{
  Isoportion lines give lines of the same portion of one of the parts,
  while isoproportion line gives lines of the same ratio between two
  parts. The isoproportion lines are Aitchison-lines in the simplex,
  while the isoportion-lines are, up to my knowledge, not.
}
\note{
Currently IsoLines only work individual plots. This is mainly due to the
fact that I have no idea, what the user interface of this function
should look like for multipanel plots. This includes philosophical
problems with the meaning of isoportions in case of marginal plots. 
  }
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{plot.acomp}}
}

\examples{
data(SimulatedAmounts)
plot(acomp(sa.lognormals))
isoPortionLines()
plot(acomp(sa.lognormals),center=TRUE)
isoPortionLines()
plot(rcomp(sa.lognormals))
isoPortionLines()
plot(acomp(sa.lognormals))
isoProportionLines()
plot(acomp(sa.lognormals),center=TRUE)
isoProportionLines()
plot(rcomp(sa.lognormals))
isoProportionLines()
}
\keyword{multivariate}
