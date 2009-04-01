extensiveCheck <- FALSE
noTestCheck  <- TRUE
if( !noTestCheck ) {
x <- runif.acomp(100,4)
y <- runif.acomp(100,4)

erg <- acompGOF.test(list(x,y))
#continue
erg
unclass(erg)
(erg <- acompGOF.test(list(x,y)))
hist(replicate(100,acompGOF.test(list(x,y))$p.value))
ks.test(dd,punif)

x <- runif.acomp(100,4)
y <- runif.acomp(100,4)
dd <- replicate(1000,acompGOF.test(list(runif.acomp(100,4),runif.acomp(100,4)))$p.value)
hist(dd)
ks.test(dd,punif)

dd <- replicate(1000,acompGOF.test(list(runif.acomp(20,4),runif.acomp(100,4)))$p.value)
hist(dd)
ks.test(dd,punif)

dd <- replicate(1000,acompGOF.test(list(runif.acomp(10,4),runif.acomp(100,4)))$p.value)
hist(dd)
ks.test(dd,punif)

dd <- replicate(1000,acompGOF.test(list(runif.acomp(10,4),runif.acomp(400,4)))$p.value)
hist(dd)
ks.test(dd,punif)

dd <- replicate(1000,acompGOF.test(list(runif.acomp(400,4),runif.acomp(10,4)))$p.value)
ks.test(dd,punif)

hist(dd)

dd <- replicate(1000,acompGOF.test(list(runif.acomp(20,4),runif.acomp(100,4)+acomp(c(1,2,3,1))))$p.value)
hist(dd)



#x <- runif.acomp(100,4)
#acompUniformityGOF.test(x)

#dd <- replicate(1000,acompUniformityGOF.test(runif.acomp(10,4))$p.value)

#hist(dd)
#shapiro.test(dd)

#x <- runif.acomp(10,4)
#y <- runif.acomp(2000,4)
#sp <- replicate(100,acompUniformityGOF.test(x,replicates=999)$p.value)
#hist(sp)

sp <- replicate(100,acompTwoSampleGOF.test(x,y,R=999)$p.value)
hist(sp)

# Gaus -Test
# Ein Stichproben
ks.test(replicate(1000,Gauss.test(rnorm(10,1,3),mean=1,sd=3)$p.value),punif)

# Zwei Stichproben
ks.test(replicate(1000,Gauss.test(rnorm(10,1,3),rnorm(23,1,3),sd=3)$p.value),punif)

# fitDirichlet
fitDirichlet(rDirichlet.acomp(1000,1:4))

# acompNormalGOF


ks.test(replicate(100,acompNormalGOF.test(rnorm.acomp(10,mean=acomp(1:4),var=ilrvar2clr(diag(1:3))))$p.value),punif)

# ccompMultinomialGOF.test

ks.test(replicate(100,ccompMultinomialGOF.test(rmultinom.ccomp(100,acomp(1:4),20))$p.value),punif)
ks.test(replicate(100,ccompMultinomialGOF.test(rpois.ccomp(100,acomp(1:4),20))$p.value),punif)

# ccompPoissonGOF.test
ks.test(replicate(100,ccompPoissonGOF.test(rmultinom.ccomp(100,acomp(1:4),20))$p.value),punif) # Ablehnung
ks.test(erg<-replicate(100,ccompPoissonGOF.test(rpois.ccomp(100,acomp(1:4),20))$p.value),punif)
erg <- ccompPoissonGOF.test(rpois.ccomp(100,acomp(1:4),20))
# PoissonGOF

ks.test(erg <- replicate(1000,PoissonGOF.test(rpois(100,20),20)$p.value),punif)
ks.test(erg <- replicate(1000,PoissonGOF.test(rpois(100,20))$p.value),punif)
erg <- replicate(1000,PoissonGOF.test(rpois(100,20),20)$p.value)
hist(erg)
(erg<- PoissonGOF.test(rpois(100,20)))

# fitSameMeanDifferentVarianceModel

x <- lapply(list(rnorm.acomp(1000,acomp(1:4),ilrvar2clr(diag(3))),
                 rnorm.acomp(2000,acomp(1:4),ilrvar2clr(diag(3))),
                 rnorm.acomp(3000,acomp(1:4),ilrvar2clr(diag(3)))),ilr)
ilr(acomp(1:4))
fitSameMeanDifferentVarianceModel(x)
ilr(acomp(1:4))



# acompNormalLocationTest.test


x<-rnorm.acomp(60,acomp(1:4),ilrvar2clr(diag(3)))
y<-rnorm.acomp(60,acomp(1:4),ilrvar2clr(diag(3)))
g<-factor(rep(c("A","B","C"),c(10,20,30)))
acompNormalLocation.test(x,g,var.equal=TRUE)
acompNormalLocation.test(x~g,list(x=x,g=g),var.equal=TRUE)
acompNormalLocation.test(split(x,g),var.equal=TRUE)

acompNormalLocation.test(x,g)
acompNormalLocation.test(x~g,list(x=x,g=g))
acompNormalLocation.test(split(x,g))

acompNormalLocation.test(x)
acompNormalLocation.test(x,y,paired=TRUE)
acompNormalLocation.test(x,y)
acompNormalLocation.test(x,y,var.equal=TRUE)

# Multi-Sample (same variance) chisq #neg 100
ks.test(erg <- replicate(100,acompNormalLocation.test(rnorm.acomp(60,acomp(1:4),ilrvar2clr(diag(3))),g,var.equal=TRUE)$p.value),punif)
hist(erg)
plot(ecdf(erg));xx<-seq(0,24,by=0.1);lines(xx,punif(xx))


# Multi-Sample (same variance) chisq #neg 1000
ks.test(erg <- replicate(1000,acompNormalLocation.test(rnorm.acomp(60,acomp(1:4),ilrvar2clr(diag(3))),g,var.equal=TRUE)$p.value),punif)
hist(erg)
plot(ecdf(erg));xx<-seq(0,24,by=0.1);lines(xx,punif(xx))

# Multi-Sample (same variance) R #pos 1000
ks.test(erg <- replicate(10000,acompNormalLocation.test(rnorm.acomp(60,acomp(1:4),ilrvar2clr(diag(3))),g,var.equal=TRUE,R=99)$p.value),punif)
hist(erg)
plot(ecdf(erg));xx<-seq(0,24,by=0.1);lines(xx,punif(xx))


#ks.test(erg <- replicate(1000,acompNormalLocation.test(rnorm.acomp(60,acomp(1:4),ilrvar2clr(diag(3))),g,var.equal=TRUE)$statistic),pchisq,df=6)
#hist(erg,freq=FALSE);xx<-seq(0,24,by=0.1);lines(xx,dchisq(xx,df=6))
#plot(ecdf(erg));xx<-seq(0,24,by=0.1);lines(xx,pchisq(xx,df=6))
 
# Multi-Sample (different variance) #neg
ks.test(erg <- replicate(1000,acompNormalLocation.test(rnorm.acomp(60,acomp(1:4),ilrvar2clr(diag(3))),g,var.equal=FALSE)$p.value),punif)
hist(erg)

# Multi-Sample (different variance) #neg
ks.test(erg <- replicate(100,acompNormalLocation.test(rnorm.acomp(60,acomp(1:4),ilrvar2clr(diag(3))),g,var.equal=FALSE,R=999)$p.value),punif)
hist(erg)


# One Sample #ok 1000 # neg 10000
ks.test(erg <- replicate(10000,acompNormalLocation.test(rnorm.acomp(60,acomp(rep(1,4)),ilrvar2clr(diag(3))))$p.value),punif)
hist(erg)

# Two Sample # neg
ks.test(erg <- replicate(100,acompNormalLocation.test(rnorm.acomp(600,acomp(rep(1,4)),ilrvar2clr(diag(3))),rnorm.acomp(600,acomp(rep(1,4)),ilrvar2clr(diag(3))))$p.value),punif)
hist(erg)

# Two Sample (equal variances) #ok 1000 # neg 10000
ks.test(erg <- replicate(10000,acompNormalLocation.test(rnorm.acomp(60,acomp(rep(1,4)),ilrvar2clr(diag(3))),rnorm.acomp(60,acomp(rep(1,4)),ilrvar2clr(diag(3))),var.equal=TRUE)$p.value),punif)
hist(erg)

# Paired # 1000 ok
ks.test(erg <-
        replicate(1000,
                  acompNormalLocation.test(rnorm.acomp(60,acomp(rep(1,4)),ilrvar2clr(diag(3))),
                                           rnorm.acomp(60,acomp(rep(1,4)),ilrvar2clr(diag(3))),
                                           paired=TRUE)$p.value),
        punif)
hist(erg)


# Internal Functions
attach(gsi)
# gsi.sortedUniforms
x<- gsi.sortedUniforms(1000)
all(x==sort(x))
ks.test(x,punif)
hist(replicate(1000,ks.test(gsi.sortedUniforms(1000),punif)$p.value))


# Likelihood Ratios

x <- ilrInv(matrix(c(
            c(1,1),
            c(-1,1),
            c(1,-1),
            c(-1,-1)
            ),byrow=TRUE,ncol=2))

acompNormalLocation.test(x,x,var.equal=FALSE)
acompNormalLocation.test(x,x,var.equal=TRUE)
acompNormalLocation.test(list(x,x,x),var.equal=TRUE)

}
