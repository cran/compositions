\name{gsi.margin}
\alias{gsi.margin}
\alias{gsi.margin.acomp}
\alias{gsi.margin.rcomp}
\alias{gsi.margin.aplus}
\alias{gsi.margin.rplus}
\title{Internal function: Compute a desired compositional margin}
\description{
  This generic function should select the selected type of margin based
  on the class of the dataset and the specified margin type. 
}
\usage{
gsi.margin(X,\dots)
gsi.margin.acomp(X,what,\dots,margin="acomp")
gsi.margin.rcomp(X,what,\dots,margin="rcomp")
gsi.margin.aplus(X,what,\dots)
gsi.margin.rplus(X,what,\dots)

}
\arguments{
  \item{X}{The dataset to take the margin from.}
  \item{what}{The indices xor column names to be kept.}
  \item{margin}{The type of marginalisation to be used. Possible values
    are: '"sub"', '"rcomp"', '"acomp"' and an index xor a name of a
    variable in the dataset.}
  \item{\dots}{ other arguments }
}
\value{
  Some marginal dataset or vector still containing the variables given
  by \code{what} and optionally one
  additional part named '"+"', '"*"' or \var{margin}.
}
\details{
}

\note{Do not use gsi.* functions directly since they are internal
  functions of the package}

\seealso{
  \code{\link{gsi}}
}
\examples{
data(SimulatedAmounts)
plot(gsi.margin(acomp(sa.lognormals5),c("Cd","Cu")))
}
\keyword{}
