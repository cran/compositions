\name{gsicall}
\alias{gsi.call}
\title{Internal functions of the compositions package}
\description{
  Calls a function with a list of arguments.
}
\usage{
gsi.call(fkt,\dots)
}

\arguments{
  \item{fkt}{The function to be called}
  \item{\dots}{The arguments to call the function with}
}
\value{
}
\details{
This is only useful in conjunction with do.call and allows to call
anonymous functions with a parameters given by a list. 
}

\note{Do not use gsi.* functions directly since they are internal
  functions of the package}

\seealso{
  \code{\link{gsi}}
}
\examples{
mypars <- list(x=3)
do.call("gsi.call",c(list(function(x){x}),mypars))
}
\keyword{}
