 \name{SimulatedAmounts}
 \docType{data}
 \alias{SimulatedAmounts}
 \alias{sa.dirichlet}
 \alias{sa.dirichlet.dil}
 \alias{sa.dirichlet.mix}
 \alias{sa.dirichlet5}
 \alias{sa.dirichlet5.dil}
 \alias{sa.dirichlet5.mix}
 \alias{sa.uniform}
 \alias{sa.uniform.dil}
 \alias{sa.uniform.mix}
 \alias{sa.uniform5}
 \alias{sa.uniform5.dil}
 \alias{sa.uniform5.mix}
 \alias{sa.lognormals}
 \alias{sa.lognormals.dil}
 \alias{sa.lognormals.mix}
 \alias{sa.lognormals5}
 \alias{sa.lognormals5.dil}
 \alias{sa.lognormals5.mix}
 \alias{sa.tnormals}
 \alias{sa.tnormals.dil}
 \alias{sa.tnormals.mix}
 \alias{sa.tnormals5}
 \alias{sa.tnormals5.dil}
 \alias{sa.tnormals5.mix}
 \alias{sa.groups}
 \alias{sa.groups.dil}
 \alias{sa.groups.mix}
 \alias{sa.groups5}
 \alias{sa.groups5.dil}
 \alias{sa.groups5.mix}
 \alias{sa.groups.area}
 \alias{sa.groups5.area}
 \alias{sa.outliers1}
 \alias{sa.outliers2}
 \alias{sa.outliers3}
 \alias{sa.outliers4}
 \alias{sa.outliers5}
 \alias{sa.outliers6}
 \alias{sa.missings}
 \alias{sa.missings5}
 \title{Simulated amount datasets}
 \description{
   Several simulated datasets intended as reference examples for various
   conceptual and statistical models of compositions and amounts.
 }
 \usage{data(SimulatedAmounts)}
 \format{Data matrices with 60 cases and 3 or 5 variables.}
 \source{The datasets are simulated for this package and are under the
   GNU Public Library Licence Version 2 or newer. }
 \details{
   The statistical analysis of amounts and compositions is set to
   discussion. Four essentially different approaches are
   provided in this package around the classes "rplus", "aplus",
   "rcomp", "acomp". There is no absolutely "right" approach, since there is
   a conection between these approaches and the processes originating the data. 
   We provide here simulated
   standard datasets and the corresponding simulation procedures
   following these several models to provide "good" analysis examples
   and to show how these models actually look like in data.
   \cr

   The data sets are simulated according to correlated lognormal
   distributions (sa.lognormals, sa.lognormal5), winsorised correlated
   normal distributions (sa.tnormals, sa.tnormal5), Dirichlet
   distribution on the simplex (sa.dirichlet, sa.dirichlet5),
   uniform distribution on the simplex  (sa.uniform, sa.uniform5), and
   a grouped dataset (sa.groups, sa.groups5) with three
   groups (given in sa.groups.area and sa.groups5.area) all distributed
   accordingly with a lognormal distribution with group-dependent means.\cr

   We can imagine that amounts evolve in nature e.g. in part of the soil they are
   diluted and transported in a transport medium, usually water, which
   comes from independent source (the rain, for instance) and this new
   composition is normalized by taking a sample of standard size.
   For each of the datasets \var{sa.X} there is a corresponding
   \var{sa.X}\code{.dil} dataset which is build by simulating exactly that process
   on the corresponding \var{sa.X} dataset . The
   amounts in the \var{sa.X}\code{.dil} are given in ppm. This idea
   of a transport medium is a major argument for a compositional
   approach, because the total amount given by the sum of the
   parts is induced by the dilution given by the medium and
   thus non-informative for the original process investigated.
   \cr

   If we imagine now these amounts flowing into a river and sedimenting, the different contributions
   are accumulated along the river and renormalized to a unit portion on
   taking samples again. For each of the dataset
   \var{sa.X}\code{.dil} there is a corresponding
   \var{sa.X}\code{.mix} dataset which is built from the corresponding
   \var{sa.X} dataset by simulating exactly that accumulation
   process. Mixing of different compositions is a major argument against
   the log based approaches (\code{\link{aplus}}, \code{\link{acomp}})
   since mixing is a highly nonlinear operation in terms of log-ratios.

 }
 \examples{
data(SimulatedAmounts)
plot.acomp(sa.lognormals)
plot.acomp(sa.lognormals.dil)
plot.acomp(sa.lognormals.mix)
plot.acomp(sa.lognormals5)
plot.acomp(sa.lognormals5.dil)
plot.acomp(sa.lognormals5.mix)

plot(acomp(sa.missings))
plot(acomp(sa.missings5))

#library(MASS)
plot.rcomp(sa.tnormals)
plot.rcomp(sa.tnormals.dil)
plot.rcomp(sa.tnormals.mix)
plot.rcomp(sa.tnormals5)
plot.rcomp(sa.tnormals5.dil)
plot.rcomp(sa.tnormals5.mix)

plot.acomp(sa.groups,col=as.numeric(sa.groups.area),pch=20)
plot.acomp(sa.groups.dil,col=as.numeric(sa.groups.area),pch=20)
plot.acomp(sa.groups.mix,col=as.numeric(sa.groups.area),pch=20)
plot.acomp(sa.groups5,col=as.numeric(sa.groups.area),pch=20)
plot.acomp(sa.groups5.dil,col=as.numeric(sa.groups.area),pch=20)
plot.acomp(sa.groups5.mix,col=as.numeric(sa.groups.area),pch=20)

plot.acomp(sa.uniform)
plot.acomp(sa.uniform.dil)
plot.acomp(sa.uniform.mix)
plot.acomp(sa.uniform5)
plot.acomp(sa.uniform5.dil)
plot.acomp(sa.uniform5.mix)

plot.acomp(sa.dirichlet)
plot.acomp(sa.dirichlet.dil)
plot.acomp(sa.dirichlet.mix)
plot.acomp(sa.dirichlet5)
plot.acomp(sa.dirichlet5.dil)
plot.acomp(sa.dirichlet5.mix)

# The data was simulated with the following commands:

#library(MASS)
dilution <- function(x) {clo(cbind(x,exp(rnorm(nrow(x),5,1))))[,1:ncol(x)]*1E6}
seqmix   <- function(x) {clo(apply(x,2,cumsum))*1E6}


vars  <- c("Cu","Zn","Pb")
vars5 <- c("Cu","Zn","Pb","Cd","Co")

sa.lognormals <- structure(exp(matrix(rnorm(3*60),ncol=3) \%*\%
                               chol(matrix(c(1,0.8,-0.2,0.8,1,
                                             -0.2,-0.2,-0.2,1),ncol=3))+
                               matrix(rep(c(1:3),each=60),ncol=3)),
                           dimnames=list(NULL,vars))

plot.acomp(sa.lognormals)
pairs(sa.lognormals)

sa.lognormals.dil <- dilution(sa.lognormals)
plot.acomp(sa.lognormals.dil)
pairs(sa.lognormals.dil)

sa.lognormals.mix <- seqmix(sa.lognormals.dil)
plot.acomp(sa.lognormals.mix)
pairs(sa.lognormals.mix)


sa.lognormals5 <- structure(exp(matrix(rnorm(5*60),ncol=5) \%*\%
                               chol(matrix(c(1,0.8,-0.2,0,0,
                                             0.8,1,-0.2,0,0,
                                             -0.2,-0.2,1,0,0,
                                             0,0,0,5,4.9,
                                             0,0,0,4.9,5),ncol=5))+
                               matrix(rep(c(1:3,-2,-2),each=60),ncol=5)),
                           dimnames=list(NULL,vars5))

plot.acomp(sa.lognormals5)
pairs(sa.lognormals5)

sa.lognormals5.dil <- dilution(sa.lognormals5)
plot.acomp(sa.lognormals5.dil)
pairs(sa.lognormals5.dil)

sa.lognormals5.mix <- seqmix(sa.lognormals5.dil)
plot.acomp(sa.lognormals5.mix)
pairs(sa.lognormals5.mix)



sa.groups.area <- factor(rep(c("Upper","Middle","Lower"),each=20))
sa.groups <- structure(exp(matrix(rnorm(3*20*3),ncol=3) \%*\%
                               chol(0.5*matrix(c(1,0.8,-0.2,0.8,1,
                                             -0.2,-0.2,-0.2,1),ncol=3))+
                               matrix(rep(c(1,2,2.5,2,2.9,5,4,2,5),
                                          each=20),ncol=3)),
                           dimnames=list(NULL,c("clay","sand","gravel")))

plot.acomp(sa.groups,col=as.numeric(sa.groups.area),pch=20)
pairs(sa.lognormals,col=as.numeric(sa.groups.area),pch=20)

sa.groups.dil <- dilution(sa.groups)
plot.acomp(sa.groups.dil,col=as.numeric(sa.groups.area),pch=20)
pairs(sa.groups.dil,col=as.numeric(sa.groups.area),pch=20)

sa.groups.mix <- seqmix(sa.groups.dil)
plot.acomp(sa.groups.mix,col=as.numeric(sa.groups.area),pch=20)
pairs(sa.groups.mix,col=as.numeric(sa.groups.area),pch=20)



sa.groups5.area <- factor(rep(c("Upper","Middle","Lower"),each=20))
sa.groups5 <- structure(exp(matrix(rnorm(5*20*3),ncol=5) \%*\%
                               chol(matrix(c(1,0.8,-0.2,0,0,
                                             0.8,1,-0.2,0,0,
                                             -0.2,-0.2,1,0,0,
                                             0,0,0,5,4.9,
                                             0,0,0,4.9,5),ncol=5))+
                               matrix(rep(c(1,2,2.5,
                                            2,2.9,5,
                                            4,2.5,0,
                                            -2,-1,-1,
                                            -1,-2,-3),
                                          each=20),ncol=5)),
                           dimnames=list(NULL,
                             vars5))

plot.acomp(sa.groups5,col=as.numeric(sa.groups5.area),pch=20)
pairs(sa.groups5,col=as.numeric(sa.groups5.area),pch=20)

sa.groups5.dil <- dilution(sa.groups5)
plot.acomp(sa.groups5.dil,col=as.numeric(sa.groups5.area),pch=20)
pairs(sa.groups5.dil,col=as.numeric(sa.groups5.area),pch=20)

sa.groups5.mix <- seqmix(sa.groups5.dil)
plot.acomp(sa.groups5.mix,col=as.numeric(sa.groups5.area),pch=20)
pairs(sa.groups5.mix,col=as.numeric(sa.groups5.area),pch=20)



sa.tnormals <- structure(pmax(matrix(rnorm(3*60),ncol=3) \%*\%
                               chol(matrix(c(1,0.8,-0.2,0.8,1,
                                             -0.2,-0.2,-0.2,1),ncol=3))+
                               matrix(rep(c(0:2),each=60),ncol=3),0),
                           dimnames=list(NULL,c("clay","sand","gravel")))

plot.rcomp(sa.tnormals)
pairs(sa.tnormals)

sa.tnormals.dil <- dilution(sa.tnormals)
plot.acomp(sa.tnormals.dil)
pairs(sa.tnormals.dil)

sa.tnormals.mix <- seqmix(sa.tnormals.dil)
plot.acomp(sa.tnormals.mix)
pairs(sa.tnormals.mix)



sa.tnormals5 <- structure(pmax(matrix(rnorm(5*60),ncol=5) \%*\%
                               chol(matrix(c(1,0.8,-0.2,0,0,
                                             0.8,1,-0.2,0,0,
                                             -0.2,-0.2,1,0,0,
                                             0,0,0,0.05,0.049,
                                             0,0,0,0.049,0.05),ncol=5))+
                               matrix(rep(c(0:2,0.1,0.1),each=60),ncol=5),0),
                           dimnames=list(NULL,
                             vars5))

plot.rcomp(sa.tnormals5)
pairs(sa.tnormals5)

sa.tnormals5.dil <- dilution(sa.tnormals5)
plot.acomp(sa.tnormals5.dil)
pairs(sa.tnormals5.dil)

sa.tnormals5.mix <- seqmix(sa.tnormals5.dil)
plot.acomp(sa.tnormals5.mix)
pairs(sa.tnormals5.mix)



sa.dirichlet <- sapply(c(clay=0.2,sand=2,gravel=3),rgamma,n=60)
colnames(sa.dirichlet) <- vars

plot.acomp(sa.dirichlet)
pairs(sa.dirichlet)

sa.dirichlet.dil <- dilution(sa.dirichlet)
plot.acomp(sa.dirichlet.dil)
pairs(sa.dirichlet.dil)

sa.dirichlet.mix <- seqmix(sa.dirichlet.dil)
plot.acomp(sa.dirichlet.mix)
pairs(sa.dirichlet.mix)



sa.dirichlet5 <- sapply(c(clay=0.2,sand=2,gravel=3,humus=0.1,plant=0.1),rgamma,n=60)
colnames(sa.dirichlet5) <- vars5

plot.acomp(sa.dirichlet5)
pairs(sa.dirichlet5)

sa.dirichlet5.dil <- dilution(sa.dirichlet5)
plot.acomp(sa.dirichlet5.dil)
pairs(sa.dirichlet5.dil)

sa.dirichlet5.mix <- seqmix(sa.dirichlet5.dil)
plot.acomp(sa.dirichlet5.mix)
pairs(sa.dirichlet5.mix)


sa.uniform   <- sapply(c(clay=1,sand=1,gravel=1),rgamma,n=60)
colnames(sa.uniform) <- vars

plot.acomp(sa.uniform)
pairs(sa.uniform)

sa.uniform.dil <- dilution(sa.uniform)
plot.acomp(sa.uniform.dil)
pairs(sa.uniform.dil)

sa.uniform.mix <- seqmix(sa.uniform.dil)
plot.acomp(sa.uniform.mix)
pairs(sa.uniform.mix)



sa.uniform5   <- sapply(c(clay=1,sand=1,gravel=1,humus=1,plant=1),rgamma,n=60)
colnames(sa.uniform5) <- vars5

plot.acomp(sa.uniform5)
pairs(sa.uniform5)

sa.uniform5.dil <- dilution(sa.uniform5)
plot.acomp(sa.uniform5.dil)
pairs(sa.uniform5.dil)

sa.uniform5.mix <- seqmix(sa.uniform5.dil)
plot.acomp(sa.uniform5.mix)
pairs(sa.uniform5.mix)

tmp<-set.seed(1400)
A <- matrix(c(0.1,0.2,0.3,0.1),nrow=2)
Mvar <- 0.1*ilrvar2clr(A \%*\% t(A))
Mcenter <- acomp(c(1,2,1))
typicalData <- rnorm.acomp(100,Mcenter,Mvar) # main population
colnames(typicalData)<-c("A","B","C")
# A dataset without outliers
sa.outliers1 <- acomp(rnorm.acomp(100,Mcenter,Mvar))
# A dataset with 10\% data with a large error in the first component
sa.outliers2 <- acomp(rbind(typicalData+rbinom(100,1,p=0.1)*rnorm(100)*acomp(c(4,1,1))))
# A dataset with a single outlier
sa.outliers3 <- acomp(rbind(typicalData,acomp(c(0.5,1.5,2))))
colnames(sa.outliers3)<-colnames(typicalData)
tmp<-set.seed(30)
rcauchy.acomp <- function (n, mean, var){
  D <- gsi.getD(mean)-1
  perturbe(ilrInv(matrix(rnorm(n*D)/rep(rnorm(n),D), ncol = D) \%*\% chol(clrvar2ilr(var))), mean)
}
# A dataset with a Cauchy type distribution
sa.outliers4 <- acomp(rcauchy.acomp(100,acomp(c(1,2,1)),Mvar/4))
colnames(sa.outliers4)<-colnames(typicalData)
# A dataset with like sa.outlier2 but a differently strong distortions
sa.outliers5 <- acomp(rbind(unclass(typicalData)+outer(rbinom(100,1,p=0.1)*runif(100),c(0.1,1,2))))
# A dataset with a second population
sa.outliers6 <- acomp(rbind(typicalData,rnorm.acomp(20,acomp(c(4,4,1)),Mvar)))

# Missings
sa.missings <- simulateMissings(sa.lognormals,dl=0.05,MAR=0.05,MNAR=0.05,SZ=0.05)
sa.missings[5,2]<-BDLvalue

sa.missings5 <- simulateMissings(sa.lognormals5,dl=0.05,MAR=0.05,MNAR=0.05,SZ=0.05)
sa.missings5[5,2]<-BDLvalue


objects(pattern="sa.*")
 }
 \references{
   %http://www.stat.boogaart.de/compositions/data \cr

   Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
     Data} Monographs on Statistics and Applied Probability. Chapman &
   Hall Ltd., London (UK). 416p.\cr

  Rehder, S. and U. Zier (2001) Letter to the Editor: Comment on 
  ''Logratio Analysis and Compositional Distance'' by J. Aitchison, C. 
  Barcel\'o -Vidal, J.A. Mart\'in-Fern\'andez and V. Pawlowsky-Glahn,
  \emph{Mathematical Geology}, \bold{33} (7), 845-848.\cr

  Zier, U. and S. Rehder (2002) Some comments on log-ratio transformation and compositional distance,
  \emph{Terra Nostra}, Schriften der Alfred Wegener-Stiftung, 03/2003\cr

 }
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\keyword{datasets}
