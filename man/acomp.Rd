\name{acomp}
\alias{acomp}
\title{Aitchison compositions}
\description{
  A class providing the means to analyse compositions in the
  philosophical framework of the Aitchison Simplex.
}
\usage{
          acomp(X,parts=1:NCOL(oneOrDataset(X)),total=1)
          }
\arguments{
  \item{X}{composition or dataset of compositions}
  \item{parts}{vector containing the indices xor names of the columns to be used}
  \item{total}{the total amount to be used, typically 1 or 100}
}
\value{
  a vector of class \code{"acomp"} representing one closed composition
  or a matrix of class \code{"acomp"} representing
  multiple closed compositions each in one row.  }

\details{
  Many multivariate datasets essentially describe amounts of D different
  parts in a whole. This has some important implications justifying to
  regard them as a scale for its own, called a
  composition. This scale  was in-depth analysed by Aitchison
  (1986) and the functions around the class \code{"acomp"} follow his
  approach.\cr
  Compositions have some important properties: Amounts are always
  positive. The amount of every part is limited to the whole. The
  absolute amount of the whole is noninformative since it is typically due
  to artifacts on the measurement procedure. Thus only relative changes
  are relevant. If the relative amount of one part
  increases, the amounts of other parts must decrease, introducing
  spurious anticorrelation (Chayes 1960), when analysed directly. Often
  parts (e.g H2O, Si) are missing in the dataset leaving the total
  amount unreported and longing for analysis procedures avoiding
  spurious effects when applied to such subcompositions. Furthermore, 
  the result of an analysis should be indepent of the units (ppm, g/l, vol.\%, mass.\%, molar
  fraction) of the dataset.
  \cr
  From these properties Aitchison showed that the
  analysis should be based on ratios or log-ratios only. He introduced 
  several transformations (e.g. \code{\link{clr}},\code{\link{alr}}),
  operations (e.g. \code{\link{perturbe}}, \code{\link{power.acomp}}),
  and a distance (\code{\link{dist}}) which are compatible
with these
  properties. Later it was found that the set of compostions equiped with
  perturbations as addition and powertransform as scalar multiplication
  and the \code{\link{dist}} as distance form a D-1 dimensional
  euclidean vector space (Billheimer, Fagan and Guttorp, 2001), which 
  can be mapped isometrically to a usual real vector space by \code{\link{ilr}} 
  (Pawlowsky-Glahn and Egozcue, 2001).
  \cr
  The general approach in analysing acomp objects is thus to performe
  classical multivariate analysis on clr/alr/ilr-transformed coordinates
  and to backtransform or display the results in such a way that they
  can be interpreted in terms of the original compositional parts.   
  \cr
  A side effect of the procedure is to force the compositions to sum up to a
  \var{total}, which is done by the closure operation \code{\link{clo}} . 
}
\seealso{
  \code{\link{clr}},\code{\link{rcomp}}, \code{\link{aplus}},
  \code{\link{princomp.acomp}}, 
  \code{\link{plot.acomp}}, \code{\link{boxplot.acomp}},
  \code{\link{barplot.acomp}}, \code{\link{mean.acomp}},
  \code{\link{var.acomp}}, \code{\link{variation.acomp}},
  \code{\link{cov.acomp}}, \code{\link{msd}}
}

\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr

  Aitchison, J, C. Barcel'o-Vidal, J.J. Egozcue, V. Pawlowsky-Glahn
  (2002) A consise guide to the algebraic geometric structure of the
  simplex, the sample space for compositional data analysis, \emph{Terra
    Nostra}, Schriften der Alfred Wegener-Stiftung, 03/2003\cr

  Billheimer, D., P. Guttorp, W.F. and Fagan (2001) Statistical interpretation of species composition,
  \emph{Journal of the American Statistical Association}, \bold{96} (456), 1205-1214\cr

  Pawlowsky-Glahn, V. and J.J. Egozcue (2001) Geometric approach to
  statistical analysis on the simplex. \emph{SERRA} \bold{15}(5), 384-398\cr

  Pawlowsky-Glahn, V. and ??? (2003) ???

  
  \url{http://ima.udg.es/Activitats/CoDaWork03}\cr
  
  \url{http://ima.udg.es/Activitats/CoDaWork05}

}
\examples{
data(SimulatedAmounts)
plot(acomp(sa.lognormals))
}
\keyword{multivariate}
