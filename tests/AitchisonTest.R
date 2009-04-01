
require(compositions)
(erg<-AitchisonDistributionIntegrals(c(-1,3,-2),ilrvar2clr(-diag(c(1,2))),grid=60))

(myvar<-with(erg, -1/2*ilrvar2clr(solve(clrvar2ilr(beta)))))
(mymean<-with(erg,myvar%*%theta))

with(erg,myvar-clrVar)
with(erg,mymean-clrMean)

res <- rAitchison(100,theta=c(1,2,3),sigma=ilrvar2clr(diag(c(1,2))))
res2<- rnorm.acomp(100,acomp(c(1,2,3)),ilrvar2clr(diag(c(1,2))))
plot(res)
plot(res2,add=TRUE,col="red")

dAitchison(res,theta=c(1,2,3),sigma=ilrvar2clr(diag(c(1,2))))

AitchisonDistributionIntegrals
