pairwisePlot <- function(X,Y,...) UseMethod("pairwisePlot",X)

                                 
pairwisePlot.default <- function(X,Y=X,...,xlab=deparse(substitute(X)),ylab=deparse(substitute(Y)),nm=c(length(Y),length(X)),panel=plot) {
  ylab
  xlab
  if( is.data.frame(X) ) {
    X<-X # data.frames are lists already
  } else if( is.list(X) ) {
    if( !missing(xlab) || is.null(names(X)) )
      names(X) <- if( length(xlab) < length(X) ) paste(xlab,1:length(X),sep="") else xlab
  } else {
    if( length(dim(X))== 0 ) {
      X <- t(t(X))
      colnames(X) <- xlab
    }
    if( is.null(colnames(X)) )
      colnames(X) <- if( length(xlab) < ncol(X) ) paste(xlab,1:ncol(X),sep="") else xlab
    cn <- 1:ncol(X)
    names(cn)<-colnames(X)
    X <- lapply(cn,function(i) X[,i])
  }
                                        #else {
                                        #stop("Unhandled type of Y in pairwisePlot:",class(Y))
                                        #}
  if( is.data.frame(Y) ){
    Y<-Y # data.frames are lists already
  } else if( is.list(Y) ) {
    if( !missing(ylab) || is.null(names(Y)) )
      names(Y) <- if( length(ylab) < length(Y) ) paste(ylab,1:length(Y),sep="") else ylab
  } else {
    if( length(dim(Y))== 0 ) {
      Y <- t(t(Y))
      colnames(Y) <- ylab
    }
    if( is.null(colnames(Y)) )
      colnames(Y) <- if( length(ylab) < ncol(Y) ) paste(ylab,1:ncol(Y),sep="") else ylab
    cn <- 1:ncol(Y)
    names(cn)<-colnames(Y)
    Y <- lapply(cn,function(i) Y[,i])
  }
                                        #else {
                                        # stop("Unhandled type of Y in pairwisePlot:",class(Y))
                                        #}
  if( !is.null(nm) ) {
    opar <- par(mfrow=nm)
    on.exit(par(opar))
  }
  for(j in 1:length(Y) )
    for(i in 1:length(X)) {
      panel(X[[i]],Y[[j]],
            xlab=names(X)[i],
            ylab=names(Y)[j],...)
    }
}
