\name{missing.compositions}
\alias{compositions.missing}
\alias{compositions.missings}
\alias{composition.missing}
\alias{composition.missings}
\alias{missings}
\alias{is.BDL}
\alias{is.SZ}
\alias{is.MAR}
\alias{is.MNAR}
\alias{is.WMNAR}
\alias{is.MNV}
\alias{is.NMV}
\alias{is.WMNAR}
\alias{is.WZERO}
\alias{has.missings}
\alias{has.missings.default}
\alias{has.missings.rmult}
\alias{SZ}
\alias{MAR}
\alias{MNAR}
\alias{BDL}
\alias{NMV}
\alias{SZvalue}
\alias{MARvalue}
\alias{MNARvalue}
\alias{BDLvalue}
\title{The policy of treatment of missing values in the "compositions" package}
\description{
  This help section discusses some general strategies of working with
  missing valuess in 
  a compositional, relative or vectorial context and shows how the various
  types of missings are represented and treated in the "compositions" package, 
  according to each strategy/class of analysis of compositions or amounts.
}\usage{
is.BDL(x,mc=attr(x,"missingClassifier"))
is.SZ(x,mc=attr(x,"missingClassifier"))
is.MAR(x,mc=attr(x,"missingClassifier"))
is.MNAR(x,mc=attr(x,"missingClassifier"))
is.NMV(x,mc=attr(x,"missingClassifier"))
is.WMNAR(x,mc=attr(x,"missingClassifier"))
is.WZERO(x,mc=attr(x,"missingClassifier"))
has.missings(x,...)
\method{has.missings}{default}(x,mc=attr(x,"missingClassifier"),...)
\method{has.missings}{rmult}(x,mc=attr(x,"missingClassifier"),...)
SZvalue
MARvalue
MNARvalue
BDLvalue
}
\arguments{
  \item{x}{A vector, matrix, acomp, rcomp, aplus, rplus object for which
    we would like to know the missing status of the entries}
  \item{mc}{A missing classifier function, giving for each value one of
    the values BDL (Below Detection Limit), SZ (Structural Zero), MAR
    (Missing at random), MNAR (Missing not at random), NMV (Not missing
    value) This functions are introduced to allow a different coding of
    the missings.}
    \item{\dots}{further generic arguments}
}
\value{A logical vector or matrix with the same shape as x stating
  wether or not the value is of the given type of missing. }
\details{
  In the context of compositional data we have to consider at least four
  types of missing and zero values:
  \item MAR (Missing at random) coded by NaN, the amount was not observed or
  is otherwise missing, in a way unrelated to its actual value. This is the "nice"
  type of missing.
  \item MNAR (Missing not at random) coded by NA, the amount was not
  observed or is otherwise missing, but it was missed in a way stochastically 
  dependent on its actual value.
  \item BDL (Below detection limit) coded by 0.0 or a negative number
  giving the detection limit; the amount was
  observed but turned out to be below the detection limit and was 
  thus rounded to zero. This is an informative version of MNAR.
  \item SZ (Structural zero) coded by -Inf, the amount is absolutely zero 
  due to structural reasons. E.g. a soil sample was dried before the analysis, 
  or the sample was preprocessed so that the fraction is removed. 
  Structural zeroes are mainly treated as MAR even though they are a kind of MNAR.
  \cr\
  Based on these basic missing types, the following extended types are defined:
  \item NMV (Not Missing Value) coded by a real number, it is just an
  actually-observed value.
  \item WMNAR (Wider MNAR) includes BDL and MNAR.
  \item WZERO (Wider Zero) includes BDL and SZ
  \cr
  Each function of type \code{is.XXX} checks the status of its argument according to 
  the XXX type of value from those above.
  \cr\cr
  Different steps of a statistical analysis and different understanding
  of the data will lead to different approaches with respect to missings and zeros.
  \cr
  In the first exploratory step, the problem is to keep the
  methods working and to make the missing structure visible in the
  analysis. The user should need as less as possible extra thinking
  about missings, an get nevertheless a true picture of the data.
  \cr
  As a second step, the analyst might want to analyse the 
  missing structure for itself. This is preliminarly provided by these functions, 
  since their result can be treated as a boolean data set in any other R function.
  \cr 
  In the later inferential steps, the problem is to get results valid
  with respect to a model. One needs to be able to look through the data
  on the true processes behind, without being distracted by artifacts
  stemming from missing values.
  \cr
  Thus, the policy on how a missing value is to be introduced into the analysis depends
  on the purpose of the analysis, the type of analysis and the model behind 
  \cr\cr
  The four philosophies work with different approaches to these problems:
  \item \code{\link{rplus}}: For positive real vectors, one can either identify BDL
  with a true 0 or impute a value relative to the detection limit, with a
  function like \code{\link{zeroreplace}}. A structural zero can either
  be seen as a true zero or as a MAR value.
  \item \code{\link{rcomp}} and \code{\link{acomp}}: For these 
  relative geometries, a true zero is an alien. Thus
  a BDL is nothing else but a small unkown value. We could either decide
  to replace the value by an imputation, or go through the whole analysis 
  keeping this lack of information in mind. 
  The main problem of imputation is that by
  closing to 1, the absolute value of the detection limit is lost, and
  the detection limit can correspond to very different portions. Raw differences 
  between \emph{all, observed or missed,} components (the ground of the rcomp geometry)
  are completely distorted by the replacement. Contrarily, log-ratios between observed
  components do not change but ratios between missed components
  dramatically depend on the replacement, e.g. typically the content of gold is some orders of
  magnitude smaller than the contend of silver even around a gold
  deposit, but far away from the deposit they both might be far under detection
  limit, leading to a ratio of 1, just because nothing was observed. SZ in compositions
  might be either seen as defining two sub-populations, one fully defined and one where 
  only a subcomposition is defined. But SZ can also
  very much be like an MAR, if only a subcomposition is measured. Thus, in general 
  we can simply understand that only a subcomposition is available, i.e. a
  projection of the true value onto a sub-space: for each observation, this sub-space 
  might be different. For MAR values, this approach
  is stricly valid, and yields unbiased estimations (because these projections are stochastically independent of the observed phenomenon). For MNAR values, the 
  projections depend on the actual value, which strictly speaking yields 
  biased estimations.
  \cr
  \code{\link{aplus}}: 
  Imputation takes place by simple replacement of the value. However
  this can lead to a dramatic change of ratios and should thus be used
  only with extra care, by the same reasons explained before. 
  \cr
  More information on how missings are actually processed can be found in the help files of each
  individual functions. 
}
\seealso{
  \code{\link{zeroreplace}}, \code{\link{rmult}}, \code{\link{ilr}},
\code{\link{mean.acomp}}, \code{\link{acomp}}, \code{\link{plot.acomp}}
}
\references{

    Boogaart, K.G. v.d., R. Tolosana-Delgado, M. Bren (2006) Concepts for
handling of zeros and missing
  values in compositional data, in E. Pirard (ed.) (2006)Proccedings of
  the IAMG'2006 Annual Conference on "Quantitative Geology from multiple
  sources", September 2006, Liege, Belgium, S07-01, 4pages,
  \url{http://www.math-inf.uni-greifswald.de/~boogaart/Publications/iamg06_s07_01.pdf}, ISBN: 978-2-9600644-0-7
  


Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr
  
  Aitchison, J, C. Barcel'o-Vidal, J.J. Egozcue, V. Pawlowsky-Glahn
  (2002) A consise guide to the algebraic geometric structure of the
  simplex, the sample space for compositional data analysis, \emph{Terra
    Nostra}, Schriften der Alfred Wegener-Stiftung, 03/2003\cr
  Billheimer, D., P. Guttorp, W.F. and Fagan (2001) Statistical interpretation of species composition,
  \emph{Journal of the American Statistical Association}, \bold{96} (456), 1205-1214\cr

  Mart\'{\i}n-Fern\'andez, J.A., C. Barcel\'o-Vidal, and V. Pawlowsky-Glahn (2003) 
  Dealing With Zeros and Missing Values in Compositional 
  Data Sets Using Nonparametric Imputation. \emph{Mathematical Geology}, \bold{35}(3)
  253-278\cr

  
 }
\examples{
require(compositions)      # load library
data(SimulatedAmounts)     # load data sa.lognormals


}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana Delgado, Matevz Bren}
\keyword{multivariate}
