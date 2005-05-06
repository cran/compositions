\name{gsi.addclass}
\alias{gsi.addclass}
\title{Internal function: give an object a derived subclass}
\description{
  This function just extends the class of an object by the given class.
}
\usage{
gsi.addclass( x , cls)
}
\arguments{
  \item{x}{The object}
  \item{cls}{the new additional class}
}
\value{
  The object x with additional class attached.
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\examples{
gsi.addclass(1:10,"goofy")
}
\keyword{}
