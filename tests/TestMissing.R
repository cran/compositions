library("compositions")
data(SimulatedAmounts)

mydata <- simulateMissings(sa.groups5,detectionlimit=0.01,knownlimit=TRUE,MAR=0.05,MNARprob=0.05,SZprob=0.05)


cdata <- acomp(mydata)
cdata
mean(cdata)
mean(acomp(sa.groups5))
plot(acomp(sa.groups5))
plot(mean(cdata),add=T,col="red")
plot(mean(acomp(sa.groups5)),add=T,col="green")
mean(cdata - mean(cdata))


cdata <- rcomp(mydata)
cdata
mean(cdata)
mean(rcomp(sa.groups5))
plot(rcomp(sa.groups5))
plot(mean(cdata),add=T,col="red",pch=20)
plot(mean(rcomp(sa.groups5)),add=T,col="green",pch=20)


cdata <- aplus(mydata)
cdata
mean(cdata)
mean(aplus(sa.groups5))
plot(aplus(sa.groups5))
plot(mean(cdata),add=T,col="red",pch=20)
plot(mean(aplus(sa.groups5)),add=T,col="green",pch=20)

cdata <- rplus(mydata)
cdata
mean(cdata)
mean(rplus(sa.groups5))
plot(rplus(sa.groups5))
plot(mean(cdata),add=T,col="red",pch=20)
plot(mean(rplus(sa.groups5)),add=T,col="green",pch=20)



