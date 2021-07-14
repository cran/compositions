\name{wrapped_functions}
\alias{anova}
\title{Standard R functions wrapped for compatibility}
\description{
  These functions of a standard R distribution (from package "base") are wrapped by naive functions
  in "compositions" with the goal to ensure their normal behavior. 
}
\usage{
anova(...)

}

\arguments{
  \item{...}{arguments passed to the original function (all!)}
}
\value{
  The same as the original function from package base (i.e. search for it with `?base::anova`)
}
\details{
  The functions documented in this page are just wrapers around base functions 
  from R that, due to a variety of reasons, need pre- or post-processing when 
  "compositions" is loaded
}

\author{Raimon Tolosana-Delgado \url{http://www.stat.boogaart.de}}
\seealso{
  \code{anova} in package "base"
}


\examples{
# anova:
data("Hydrochem") # load data
Z = acomp(Hydrochem[,7:19]) # select composition
Hydrochem$compo = Z # attach to dataset
md = lm(alr(compo)~log(H), data=Hydrochem) # fit model
anova(md)   # anova test
}
\keyword{multivariate}
