\name{plot.acomp}
\alias{plot.acomp}
\alias{plot.rcomp}
\title{Displaying compositions in ternary diagrams}
\description{
}
\usage{
\method{plot}{acomp}(x,...,labels=colnames(X),cn=colnames(X),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),margin="acomp",add=FALSE,triangle=!add,col=par("col"))
\method{plot}{rcomp}(x,...,labels=colnames(X),cn=colnames(X),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),margin="rcomp",add=FALSE,col=par("col"))
          }
\arguments{
  \item{x}{a dataset of a compositional class}
  \item{\dots}{further graphical parameters passed (see
    \code{\link[graphics]{par}})}
  \item{margin}{The type of marginalisation to be computed, when
    displaying the individual panels. Possible values are: \code{"acomp"},
    \code{"rcomp"} and any of the variable names/column numbers in the
    composition. If one of the columns is selected each panel displays a
    subcomposition given by the row part, the column part and
    the given part. If one of the classes is given the corresponding
    margin \code{\link{acompmargin}}  or  \code{\link{rcompmargin}} is
    used. }
  \item{add}{a logical indicating whether the information should just
    be added to an existing plot. In case of false a new plot is
    created.}
  \item{triangle}{A logical indicating whether the triangle should be
    drawn. }
  \item{col}{The color to plot the data.}
  \item{labels}{The names of the parts}
  \item{cn}{The names of the parts to be used in a single
    panel. Internal use only.}
  \item{aspanel}{Logical indicating that only a single panel should be
    drawn and not the whole plot. Internal use only.}
  \item{id}{A logical. If true one can identify the points like with the
    \code{\link[graphics]{identify}} command. }
  \item{idlabs}{A character vector providing the labels to be used with
    the identification, when \code{id=TRUE}}
  \item{idcol}{color of the idlabs-labels}
  \item{center}{a logical indicating whether a the data should be
    centered prior to the plot. Centering is done in the choosen
    philosophy. See \code{\link{scale}}}
  \item{scale}{a logical indicating whether a the data should be
    scaled prior to the plot. Scaling is done in the choosen
    philosophy. See \code{\link{scale}}}
  \item{pca}{A logical indicating whether the first principle component
    should be displayed in the plot. Currently direction of the
    principle component
    of the displayed subcomposition is displayed as a line. Later on a
    the principle componenent of the whole dataset should be displayed.}
  \item{col.pca}{The color to draw the principle component.}
}
\details{
  The data is displayed in ternary diagrams. This does not work for
  two part compositions. Compositions of three parts are displayed
  in a single ternary diagram. For compositions of more than three
  components, the data is arrange in a scatterplot matrix through the
  command \code{\link[graphics]{pairs}}.
  \cr
  The third component in each
  of the panels is than choosen according to setting of
  \code{margin=}. Possible values of \code{margin=} are:
  \code{"acomp"},
  \code{"rcomp"} and any of the variable names/column numbers in the
  composition. If one of the columns is selected each panel displays a
  subcomposition given by the row part, the column part and
  the given part. If one of the classes is given the corresponding
  margin \code{\link{acompmargin}}  or  \code{\link{rcompmargin}} is
  used.
  \cr
  Ternary diagrams can be read in multiple ways. Each corner of the
  triangle corresponds to a composition only containing the single part
  displayed in that corner. Points on the edges correspond to
  compositions only containing the parts in the adjacent corners. The
  relative amounts are displayed by the distance to the opposite
  corner. The individual portions of general points can be infered by
  imaginatorily drawing a line parallel to the edge opposite to the
  corner of the part of interest through the point. The portion of the
  part of intrest is constant along the line. Thus we can read it 
  on both crossing points of the line with the edges.
  \cr
  Relative portions of two parts can be inferred by imaginatorily
  drawing a line through the point and the corner of the unimportant
  component. This line intersects the edge between the two components
  of interest in the composition with the same relative portion of the
  two remaining components.
  \cr
  Exactly the lines parallel to one of the edges or going through one of
  the corners are straight lines as well in Aitchison and as in real
  geometry. They remain straight under an arbitrary perturbation. 

} \references{
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
  
  \url{http://ima.udg.es/Activitats/CoDaWork03}\cr
  
  \url{http://ima.udg.es/Activitats/CoDaWork05}

}

\seealso{
  \code{\link{plot.aplus}},
  \code{\link{qqnorm.acomp}},\code{\link{boxplot.acomp}} }

\examples{
data(SimulatedAmounts)
plot(acomp(sa.lognormals))
plot(rcomp(sa.lognormals))
plot(acomp(sa.lognormals5),pca=TRUE)
plot(rcomp(sa.lognormals5),pca=TRUE)
}
\keyword{multivariate}
