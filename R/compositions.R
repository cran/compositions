# (C) 2005 by Gerald van den Boogaart, Greifswald
# License: GPL version 2 or newer

gsi.plain <- function(x) {
  if( is.data.frame(x) )
    unclass(data.matrix(x))
  else 
    unclass(x)
}

gsi.simshape <- function(x,oldx) {
  if(length(dim(oldx))>=2 )
    oneOrDataset(x)
  else if( length(dim(oldx)) == 1 )
    structure(c(x),dim=length(x))
  else 
    c(drop(x))
}

gsi.diagExtract <- function(x) {
  if( length(x) > 1 )
    diag(x)
  else
    c(x)
}

gsi.diagGenerate <- function(x) {
  if( length(x) > 1 )
    diag(x)
  else
    matrix(x)
}

gsi.getD  <- function(x) ncol(oneOrDataset(x))
gsi.getN  <- function(x) nrow(oneOrDataset(x))   

# NA   = MAR  =  (Missing at random) = not reported measured
# -Inf = MNAR =  (Missing not at random) = danger
# NaN  = SZ   =  (Structural Zero) = It makes no sense to speak about
# 0.0    = BDL  =  (Below detection limit) = Observed Zero


MARvalue <- NaN
MNARvalue<- NA
SZvalue  <- -Inf
BDLvalue <- 0.0   # Use with care, negative values specifiy the detection limit

gsi.recodeM2C <- function(x,y=x,BDL,SZ,MAR,MNAR,NMV) {
  if( !missing(BDL) ) y[is.BDL(x)]<-BDL
  if( !missing(SZ)  ) y[is.SZ(x)] <-SZ
  if( !missing(MAR) ) y[is.MAR(x)]<-MAR
  if( !missing(MNAR)) y[is.MNAR(x)]<-MNAR
  if( !missing(NMV) ) y[is.NMV(x)]<-NMV
  y
}

#gsi.recodeM2Clean <- function(x,y=x,BDL=NaN,SZ=NaN,MAR=NaN,MNAR=NA) {
#  y[is.BDL(x)]<-BDL
#  y[is.SZ(x)] <-SZ
#  y[is.MAR(x)]<-MAR
#  y[is.MNAR(x)]<-MNAR
#  y[is.NMV(x)]<-NMV
#  y
#}



gsi.recodeC2M <- function(x,y=x,na,nan,ninf,inf,neg,zero,pos) {
  if( !missing(na)  ) y[is.na(x)&!is.nan(x)]<-na
  if( !missing(nan) ) y[is.nan(x)]          <-nan
  if( !missing(ninf)) y[is.infinite(x)&x<0] <-ninf
  if( !missing(inf) ) y[is.infinite(x)&x>0] <-inf
  if( !missing(neg) ) y[is.finite(x)&x<0]   <-neg
  if( !missing(zero)) y[is.finite(x)&x==0]  <-zero
  if( !missing(pos))  y[is.finite(x)&x>0]   <-pos
  y
}

gsi.recodeM2Clean <- function(x,y=x,BDL=NaN,SZ=NaN,MAR=NaN,MNAR=NA,NMV) {
  if( !missing(BDL) ) y[is.BDL(x)]<-BDL
  if( !missing(SZ)  ) y[is.SZ(x)] <-SZ
  if( !missing(MAR) ) y[is.MAR(x)]<-MAR
  if( !missing(MNAR)) y[is.MNAR(x)]<-MNAR
  if( !missing(NMV) ) y[is.NMV(x)]<-NMV
  y
}


gsi.cleanR <- function(x) {
  x[!is.finite(x)]<-0.0
  x
}

gsi.eq <-  function(x,y) {
  if( is.null(y) ) return(is.null(x)) # null
  if( is.finite(y) ) {           
    if( is.infinite(1/y) & (1/y)<0 )  # -0
      return(is.infinite(1/x) & (1/x)<0)
      else
        return(is.finite(x) & x==y)     # Zahlencodes 
  }
  if( is.nan(y) ) return(is.nan(x))   # NaN
  if( is.infinite(y)&y>0) return(is.infinite(x)&x>0) # Inf
  if( is.infinite(y)&y<0) return(is.infinite(x)&x<0) # -Inf
  if( is.na(y) ) return(is.na(y))                    # NA
  stop("Unkown comparison type ",y)                  # Was wurde vergessen?
}


names.acomp <- function(x) colnames(oneOrDataset(x))
names.rcomp <- names.acomp
names.aplus <- names.acomp
names.rplus <- names.acomp
names.rmult <- names.acomp

"names<-.acomp" <- "names<-.rcomp" <- "names<-.aplus" <- "names<-.rplus" <- "names<-.rmult" <-
  function(x,value) {
    if(is.matrix(x)) {
      colnames(x) <- value
      x
    }
    else
      NextMethod("names",x,value=value)
  }


groupparts <- function(x,...) UseMethod("groupparts",x)

groupparts.rcomp <- function(x,...,groups=list(...)) {
                                        # BDL=SZ->0, MAR->MAR, MNAR->MNAR
  x <- rmult(gsi.recodeM2C(x,clo(x),BDL=0.0,SZ=0.0,MAR=NaN,MNAR=NA))
  usedparts <- unique(unlist(lapply(groups,function(i) {
    if( is.character(i) ) {
      parts <- match(i,names(x))
      if( any(is.na(parts)))
        stop("Unknown part",i[is.na(parts)])
      parts
    } else i
  })))
  otherparts <- (1:gsi.getD(x))[-usedparts]
  if( length(otherparts) >0 ) {
    names(otherparts) <- names(x)[otherparts]
    groups <- c(groups,otherparts)
  }
  erg <- sapply(groups,function(idx) {
    ss <- rplus(x,idx)
    ss %*% rep(1,gsi.getD(ss))
  })
  rcomp(gsi.recodeC2M(erg,na=MNARvalue,nan=MARvalue))
}

groupparts.rplus <- function(x,...,groups=list(...)) {
                                        # BDL=SZ->0, MAR->MAR, MNAR->MNAR
  x <- rmult(gsi.recodeM2C(x,BDL=0.0,SZ=0.0,MAR=NaN,MNAR=NA))
  usedparts <- unique(unlist(lapply(groups,function(i) {
    if( is.character(i) ) {
      parts <- match(i,names(x))
      if( any(is.na(parts)))
        stop("Unknown part",i[is.na(parts)])
      parts
    } else i
  })))
  otherparts <- (1:gsi.getD(x))[-usedparts]
  if( length(otherparts) >0 ) {
    names(otherparts) <- names(x)[otherparts]
    groups <- c(groups,otherparts)
  }
  erg <- sapply(groups,function(idx) {
    ss <- rplus(x,idx)
    ss %*% rep(1,gsi.getD(ss))
  })
  rplus(gsi.recodeC2M(erg,na=MNARvalue,nan=MARvalue))
}

groupparts.acomp <- function(x,...,groups=list(...)) {
                                        # BDL: BDL, NA: NA, 0: 0
  x <- rmult(gsi.recodeM2C(x,clo(x),BDL=NaN,SZ=NA,MAR=Inf,MNAR=NaN))
  #SZ <- is.SZ(x)               # keep regardless of the rest NA
  #MNAR <- is.MNAR(x)|is.BDL(x) # keep if no SZ               NaN
  #MAR  <- is.MAR(x)            # keep if no SZ or MNAR are in the way Inf
  usedparts <- unique(unlist(lapply(groups,function(i) {
    if( is.character(i) ) {
      parts <- match(i,names(x))
      if( any(is.na(parts)))
        stop("Unknown part",i[is.na(parts)])
      parts
    } else i
  })))
  otherparts <- (1:gsi.getD(x))[-usedparts]
  if( length(otherparts) >0 ) {
    names(otherparts) <- names(x)[otherparts]
    groups <- c(groups,otherparts)
  }
  erg <- sapply(groups,function(idx) {
    ss <- aplus(x,idx)
    if( is.matrix(ss) )
      gsi.geometricmean.row(ss)
    else
      gsi.geometricmean(ss)
  })
  acomp(gsi.recodeC2M(erg,na=SZvalue,nan=MNARvalue,inf=MARvalue))
  
#  SZ   <- is.na(x)&!is.na(x)   # keep regardless of the rest NA
#  MNAR <- is.nan(x)            # keep if no SZ               NaN
#  MAR  <- !is.finite(x)&!is.na(x) # keep if no SZ or MNAR are in the way Inf
}

groupparts.aplus <- function(x,...,groups=list(...)) {
  x <- rmult(gsi.recodeM2C(x,clo(x),BDL=NaN,SZ=NA,MAR=Inf,MNAR=NaN))
  usedparts <- unique(unlist(lapply(groups,function(i) {
    if( is.character(i) ) {
      parts <- match(i,names(x))
      if( any(is.na(parts)))
        stop("Unknown part",i[is.na(parts)])
      parts
    } else i
  })))
  otherparts <- (1:gsi.getD(x))[-usedparts]
  if( length(otherparts) >0 ) {
    names(otherparts) <- names(x)[otherparts]
    groups <- c(groups,otherparts)
  }
  erg <- sapply(groups,function(idx) {
    ss <- aplus(x,idx)
    if( is.matrix(ss) )
      gsi.geometricmean.row(ss)
    else
      gsi.geometricmean(ss)
  })
  aplus(gsi.recodeC2M(erg,na=SZvalue,nan=MNARvalue,inf=MARvalue))
}

# groupparts(x,G1=c("Cd","S"),G2=c("Co","Ni"),G3=c("As","F"))

clo <- function(X,parts=1:NCOL(oneOrDataset(X)),total=1,
                detectionlimit=attr(X,"detectionlimit"),
                BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL,
                storelimit=!is.null(attr(X,"detectionlimit"))) {
  X <- gsi.plain(X)
  # Collect the parts
  parts  <- unique(parts)
  if( is.character(parts) ) {
    partsn <- match(parts,colnames(X))
    if( any(is.na(partsn)) )
      stop("Unknown variable name",parts[is.na(partsn)])
    parts <- partsn
  }
  nparts <- length(parts)
  Xn <- gsi.plain(oneOrDataset(X))[,parts,drop=FALSE]
  drop <- length(dim(X)) < 2
  #if( any(na.omit(c(Xn)<0)) )
  #  stop("Negative values are not valid for amounts")
  # Processing of missings
  iMAR <- if( !is.null(MAR) ) gsi.eq(Xn,MAR) else FALSE
  iMNAR<- if( !is.null(MNAR) ) gsi.eq(Xn,MNAR) else FALSE
  iSZ  <- if( !is.null(SZ) ) gsi.eq(Xn,SZ) else FALSE
  iBDL <- if( !is.null(BDL) ) gsi.eq(Xn,BDL) else FALSE
  if( is.null(detectionlimit) ) {
    if( any(iBDL) )
      Xn[iBDL]<- BDLvalue
  } else if( is.matrix(detectionlimit) ) {
    if( nrow(Xn)!=nrow(detectionlimit) | ncol(Xn)!=ncol(detectionlimit) )
      warning("Matrix of Detectionlimits does not fit x")
    Xn <- ifelse( is.finite(detectionlimit) & detectionlimit>=0,
                ifelse( (is.finite(Xn) & X <= detectionlimit)|iBDL ,
                       -detectionlimit, Xn ),
                Xn)
  } else if( length( detectionlimit) > 1 ) {
    if( ncol(Xn)!=length(detectionlimit)  )
      warning("Length of Detectionlimits does not fit x")
    detectionlimit <- outer(rep(1,nrow(Xn),detectionlimit))
    Xn <- ifelse( is.finite(detectionlimit) & detectionlimit>=0,
                ifelse( (is.finite(Xn) & X <= detectionlimit)|iBDL ,
                       -detectionlimit, Xn ),
                Xn)
  } else if( is.finite(detectionlimit) && detectionlimit > 0 ) {
    Xn <- ifelse((Xn<=detectionlimit&Xn>=0)|iBDL,-detectionlimit,Xn)
  } else
    Xn <- ifelse(iBDL,BDLvalue,Xn)
  if( any(iMAR) ) Xn[iMAR] <- MARvalue
  if( any(iMNAR)) Xn[iMNAR]<- MNARvalue
  if( any(iSZ) )  Xn[iSZ]  <- SZvalue
  # Make the sum 1 ignoring missings
  scaling <- 1
  if( length(total) > 1 || !is.na(total) ) {
    nas <- !is.NMV(Xn)&!is.BDL(Xn) # Missings are not included in closing
    bdl <- is.BDL(Xn)              # BDLs are accordingly scaled 
    naValues <- Xn[nas]
    Xn[nas]<-0
    s <- c(ifelse(bdl,0,Xn) %*% rep(1,nparts))
    scaling <-  matrix(rep(s/total,nparts),ncol=nparts)
    Xn  <- Xn / scaling
    Xn[nas] <- naValues
  }
  erg <- gsi.simshape(Xn,X)
  if( storelimit) {
    if( length(detectionlimit) == 1 )
      detectionlimit <- matrix(detectionlimit,nrow=nrow(Xn),ncol=ncol(Xn))
    detectionlimit/scaling
    detectionlimit <- gsi.simshape(detectionlimit[,parts,drop=FALSE],X)
    attr(erg,"detectionlimit") <- detectionlimit
  }
  erg
}




acomp <- function(X,parts=1:NCOL(oneOrDataset(X)),total=1,warn.na=FALSE,detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL) {

  X <-  structure(clo(X,parts,total),class="acomp")
  if( any(is.NMV(X)&X<0)) 
    warning("Composition has negative values")
  if( warn.na ) {
    if( any(is.SZ(X))) 
      warning("Composition has structural zeros")
    if( any(is.MAR(X) | is.MNAR(X)))
      warning("Composition has missings")
    if( any(is.BDL(X)) )
      warning("Composition has values below detection limit")
  }
  X
}

rcomp <- function(X,parts=1:NCOL(oneOrDataset(X)),total=1,warn.na=FALSE,detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL) {
  X <-  structure(clo(X,parts,total,detectionlimit=detectionlimit,BDL=BDL,MAR=MAR,MNAR=MNAR,SZ=SZ),class="rcomp")
  if( warn.na ) {
    if( any(is.SZ(X))) 
      warning("Composition has structural zeros")
    if( any(is.MAR(X) | is.MNAR(X)))
      warning("Composition has missings")
    #if( any(is.BDL(X)) )
     # warning("Composition has values below detection limit")
  }
  X
}


aplus <- function(X,parts=1:NCOL(oneOrDataset(X)),total=NA,warn.na=FALSE,detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL) {
  X <- gsi.simshape(clo(X,parts,total,detectionlimit=detectionlimit,BDL=BDL,MAR=MAR,MNAR=MNAR,SZ=SZ),X)
  if( warn.na ) {
    if( any(is.SZ(X))) 
      warning("aplus has structural zeros")
    if( any(is.MAR(X) | is.MNAR(X)))
      warning("aplus has missings")
    if( any(is.BDL(X)) )
      warning("aplus has values below detection limit")
  }
  class(X) <-"aplus"
  X
}

rplus <- function(X,parts=1:NCOL(oneOrDataset(X)),total=NA,warn.na=FALSE,detectionlimit=NULL,BDL=NULL,MAR=NULL,MNAR=NULL,SZ=NULL) {
  X <- gsi.simshape(clo(X,parts,total,detectionlimit=detectionlimit,BDL=BDL,MAR=MAR,MNAR=MNAR,SZ=SZ),X)
  if( warn.na ) {
    if( any(na.omit(c(X)==0)) )
      warning("rplus has structural zeros")
    if( any(is.na(c(X)) & ! is.nan(c(X))))
      warning("rplus has missings")
    if( any(is.nan(c(X))))
      warning("rplus has values below detection limit")
  }
  class(X) <-"rplus"
  X
}

rmult <- function(X,parts=1:NCOL(oneOrDataset(X)),
                  orig=attr(X,"orig"),missingProjector=attr(X,"missingProjector")) {
  X <- gsi.simshape(oneOrDataset(X)[,parts,drop=FALSE],X)
  attr(X,"orig") <- orig
  attr(X,"missingProjector")<-missingProjector
  class(X) <-"rmult"
  X
}

print.rmult <- function(x,...) {
  Odata <- attr(x,"orig")
  if( ! is.null(Odata) )
    attr(x,"orig") <- missingSummary(Odata)
  mp <- attr(x,"missingProjector")
  if( ! is.null(mp) )
    attr(x,"missingProjector") <- dim(mp)
  NextMethod(x,...)
}
gsi2.invperm <- function(i,n){
  i <- unique(c(i,1:n))
  j <- numeric(length(i))
  j[i]<-1:length(i)
  j
}


rcompmargin <- function(X,d=c(1,2),name="+",pos=length(d)+1) {
  X <- rcomp(X)
  drop <- length(dim(X)) < 2
  if( mode(d)=="character" )
    d <- match(d,colnames(X))
  X <- oneOrDataset(gsi.plain(X))
  d <- unique(d)
  if( NCOL(X) <= length(d) )
    return(rcomp(X))
  else if( NCOL(X) == length(d) +1)
    return( rcomp(cbind(X[,d,drop=FALSE],X[,-d,drop=FALSE]) ))
  Xm <- gsi.recodeM2C(X[,-d,drop=FALSE],BDL=0.0,SZ=0.0,MAR=NaN,MNAR=NA)
  rest =gsi.recodeC2M(Xm %*% rep(1,NCOL(Xm)),zero=BDLvalue,nan=MARvalue,na=MNARvalue)
  tmp <- rcomp(cbind(rest=rest ,X[,d,drop=FALSE] ))
  if( !is.null(colnames(tmp)) )
    colnames(tmp)[1]<-name
  if( pos != 1 )
    tmp <- tmp[,gsi2.invperm(pos,ncol(tmp))]
  if( drop )
    tmp <- drop(tmp)
  rcomp(tmp)
}

acompmargin <- function(X,d=c(1,2),name="*",pos=length(d)+1) {
  drop <- length(dim(X)) < 2
  if( mode(d)=="character" )
    d <- match(d,colnames(X))
  X <- oneOrDataset(gsi.plain(X))
  d <- unique(d)
  if( NCOL(X) <= length(d) )
    return(X)
  else if( NCOL(X) == length(d) +1)
    return( cbind(X[,d,drop=FALSE],X[,-d,drop=FALSE]) )
  Xm <- X[,-d,drop=FALSE]
  Xm <- gsi.recodeM2C(Xm,log(Xm),BDL=-Inf,MAR=NaN,MNAR=NA)
  rest <- Xm %*% rep(1/NCOL(Xm),NCOL(Xm))
  rest <- gsi.recodeC2M(rest,exp(rest),inf=BDLvalue,nan=MARvalue,na=MNARvalue)
  tmp <- acomp(cbind(rest=rest ,X[,d,drop=FALSE] ))
  if( !is.null(colnames(tmp)) )
    colnames(tmp)[1]<-name
  if( pos != 1 )
    tmp <- tmp[,gsi2.invperm(pos,ncol(tmp))]
  if( drop )
    tmp <- drop(tmp)
  acomp(tmp)
}


oneOrDataset <- function(W,B=NULL) {
  W <- gsi.plain(W)
  if( missing(B) || length(dim(B))!= 2 ) {
    if( length(dim(W)) == 2) {
      return( W )
    }
    else {
      tmp <- matrix(c(W),nrow=1)
      colnames(tmp) <- names(W)
      return(tmp)
    }
  } else {
    if( length(dim(W)) == 2) {
      return( W )
    }
    else {
      tmp <- matrix(c(W),nrow=NROW(B),ncol=length(W),byrow=TRUE)
      colnames(tmp)<- names(W)
      return(tmp)
    }
  }
}

gsi.geometricmean <- function(x,...) {
    exp(mean(log(c(unclass(x))),...))
}

gsi.geometricmean.row <- function(x,...) apply(x,1,gsi.geometricmean,...)
gsi.geometricmean.col <- function(x,...) apply(x,2,gsi.geometricmean,...)


geometricmean <- function(x,...) {
  if( any(na.omit(x==0)) )
    0
  else
    exp(mean(log(unclass(x)[is.finite(x)&x>0]),...))
}

geometricmean.row <- function(x,...) apply(x,1,geometricmean,...)
geometricmean.col <- function(x,...) apply(x,2,geometricmean,...)

mean.col <- function( x , ... , na.action=get(getOption("na.action"))) {
  apply(oneOrDataset(x),2,function(x,...) mean(x,na.action=na.action),...)
}

mean.row <- function( x , ... , na.action=get(getOption("na.action"))) {
  apply(oneOrDataset(x),1,function(x,...) mean(x,na.action=na.action),...)
}


totals <- function( x , ... ) UseMethod("totals",x)

totals.acomp <- function(x,...,missing.ok=TRUE) {
  x <- oneOrDataset(x)
  if( missing.ok )
    x <- gsi.recodeM2C(x,BDL=0.0,SZ=0.0,MAR=0.0,MNAR=0.0)
  else 
    x <- gsi.recodeM2C(x,BDL=0.0,SZ=0.0,MAR=NaN,MNAR=NA)
  erg <- gsi.recodeC2M(apply(x,1,sum,...),zero=BDLvalue,nan=MARvalue,na=MNARvalue)
}

totals.aplus <- totals.acomp
totals.rcomp <- totals.acomp
totals.rplus <- totals.acomp

gsi.svdsolve <- function(a,b,...,cond=1E-10) {
  s <- svd(a,...)
  lambda1 <- s$d[1] 
  drop(s$v %*% (gsi.diagGenerate(ifelse(s$d<lambda1*cond,0,1/s$d)) %*% (t(s$u) %*% b)))
  
}

missingProjector <- function(x,...,by="s") UseMethod("missingProjector")

missingProjector.acomp <- function(x,has=is.NMV(x),...,by="s") {
  if( is.tensor(has) ) {
  } else if( is.matrix(has) ) {
    has <- as.tensor(has)
    names(has) <- c(by,"i")
  } else
    has <- to.tensor(c(has),c(i=length(has)))
  iii <- function(x) {a<-ifelse(unclass(x)==0,0,1/unclass(x));attributes(a)<-attributes(x);a}
  reorder.tensor( diag.tensor(has,by=by,mark="'") - 
                 einstein.tensor(has,mark(has,mark="'",by=by),
                                 diag=iii(margin.tensor(has,by=by)),
                                 by=by))
  #  gsi.diagGenerate(as.numeric(has))-(has %o% has)/sum(has) 
}

missingProjector.rcomp <- function(x,has=!(is.MAR(x)|is.MNAR(x)),...,by="s") {
  if( is.tensor(has) ) {
  } else  if( is.matrix(has) ) {
    has <- as.tensor(has)
    names(has) <- c(by,"i")
  } else
    has <- to.tensor(c(has),c(i=length(has)))
  iii <- function(x) {a<-ifelse(unclass(x)==0,0,1/unclass(x));attributes(a)<-attributes(x);a}
  reorder.tensor( diag.tensor(has,by=by,mark="'") -
                 einstein.tensor(has,has[[i=~"i'"]],
                                 diag=iii(margin.tensor(has,by=by)),
                                 by=by))
                 
  #  gsi.diagGenerate(as.numeric(has))-(has %o% has)/sum(has) 
}


missingProjector.aplus <- function(x,has=is.NMV(x),...,by="s") {
  if( is.tensor(has) ) {
  } else if( is.matrix(has) ) {
    has <- as.tensor(has)
    names(has) <- c(by,"i")
  } else
    has <- to.tensor(c(has),c(i=length(has)))
  reorder.tensor( diag.tensor(has,by=by,mark="'"), by=by) 
  #  gsi.diagGenerate(as.numeric(has)) 
}

missingProjector.rplus <- function(x,has=!(is.MAR(x)|is.MNAR(x)),...,by="s") {
  if( is.tensor(has) ) {
  } else if( is.matrix(has) ) {
    has <- as.tensor(has)
    names(has) <- c(by,"i")
  } else
    has <- to.tensor(c(has),c(i=length(has)))
  reorder.tensor( diag.tensor(has,by=by,mark="'"), by=by) 
  #  gsi.diagGenerate(as.numeric(has)) 
}
  
missingProjector.rmult <- function(x,has=is.finite(x),...,by="s") {
  if( !is.null(attr(x,"missingProjector") ) )
     return(attr(x,"missingProjector"))
  if( !is.null(attr(x,"orig") ) )
     return(missingProjector(attr(x,"orig") ))
  if( is.tensor(has) ) {
  } else if( is.matrix(has) ) {
    has <- as.tensor(has)
    names(has) <- c(by,"i")
  } else
    has <- to.tensor(c(has),c(i=length(has)))
  reorder.tensor( diag.tensor(has,by=by,mark="'"), by=by) 
  #  gsi.diagGenerate(as.numeric(has)) 
}


sumMissingProjector <- function(x,...) UseMethod("sumMissingProjector")

sumMissingProjector.acomp <- function(x,has=is.NMV(x),...) {
  n <- nrow(has)
  m <- ncol(has)
  canz <- rep(1,n) %*% has 
  ranz <- c(has %*% rep(1,m))
  erg <- diag(c(canz)) -  t(has/ranz) %*% has 
  erg
}

sumMissingProjector.aplus <- function(x,has=is.NMV(x),...) {
  n <- nrow(has)
  canz <- rep(1,n) %*% has 
  erg <- diag(c(canz))
  erg
}

sumMissingProjector.rcomp <- function(x,has=!(is.MAR(x)|is.MNAR(x)),...) {
  n <- nrow(has)
  m <- ncol(has)
  canz <- rep(1,n) %*% has 
  ranz <- c(has %*% rep(1,m))
  erg <- diag(c(canz)) -  t( has / ranz ) %*% has 
  erg
}

sumMissingProjector.rplus <- function(x,has=!(is.MAR(x)|is.MNAR(x)),...) {
  n <- nrow(has)
  canz <- rep(1,n) %*% has 
  erg <- diag(c(canz))
  erg
}

sumMissingProjector.rmult <- function(x,has=is.finite(x),...) {
  if( !is.null(mp<- attr(x,"missingProjector"))  )
    return( mp %e% one.tensor(c(s=dim(mp)["s"])))
  if( !is.null(orig<-attr(x,"orig")))
     return(sumMissingProjector(orig))
  n <- nrow(has)
  canz <- rep(1,n) %*% has 
  erg <- diag(canz)
  erg
}

  
gsi.csum <- function(x) {c(rep(1,nrow(x)) %*% ifelse(is.finite(x),x,0))}
gsi.rsum <- function(x) {c(ifelse(is.finite(x),x,0)%*%rep(1,ncol(x))) }


mean.aplus<- mean.rplus <- mean.rcomp <- mean.acomp <- function( x, ... ) {
  cdt.inv(mean(cdt(x),...),x)
}
 
#mean.acomp <- function( x, ... ) {
#  if( has.missings(x) ) 
#    clr.inv(gsi.svdsolve(sumMissingProjector(x),gsi.csum(gsi.cleanR(clr(x)))))
#  else
#    clr.inv(mean(clr(x),...))
#}


#mean.rcomp <- function( x,... ) {
#  if( has.missings(x) ) 
#    cpt.inv(gsi.svdsolve(sumMissingProjector(x), gsi.csum(gsi.cleanR(cpt(x))) ) )
#  else
#    cpt.inv(mean.col(cpt(x),...))
#}

#mean.aplus <- function( x,... ) {
#  if( has.missings(x) ) 
#    ilt.inv(gsi.svdsolve(sumMissingProjector(x), gsi.csum(gsi.cleanR(ilt(x))) ) )
#  else
#    ilt.inv(mean.col(ilt(x),...))
#}

#mean.rplus <- function( x,...) {
#  if( has.missings(x) ) 
#    iit.inv(gsi.svdsolve(sumMissingProjector(x), gsi.csum(gsi.cleanR(iit(x))) ) )
 # else
#    iit.inv(mean.col(iit(x),...))
#}

mean.rmult <- function( x,...,na.action=get(getOption("na.action"))) {
  if( has.missings(x) ) 
    rmult(gsi.svdsolve(sumMissingProjector(x), gsi.csum(x,...)))
  else
    rmult(mean.col(unclass(x),...,na.action=na.action))
}


clr2ilr <- function( x , V=ilrBase(x) ) {
  gsi.simshape(gsi.recodeC2M(oneOrDataset(x),ninf=0,nan=0,na=0,inf=0) %*% V , x)
}

ilr2clr <- function( z , V=ilrBase(z=z) ) {
  gsi.simshape( oneOrDataset(z) %*% t(V) , z)
}


clrvar2ilr <- function( varx , V=ilrBase(D=ncol(varx)) ) {
  t(V) %*% varx %*% V
}

ilrvar2clr <- function( varz , V=ilrBase(D=ncol(varz)+1) ) {
  V %*% varz %*% t(V)
}



var         <- function(x,...) UseMethod("var",x)
var.default <- stats::var
formals(var.default) <- c(formals(var.default),alist(...= ))


  
var.acomp   <- function(x,y=NULL,...) {
  if( has.missings(x) ||  has.missings(y) ) 
    return(var(cdt(x),cdt(y),...))
  else
    return(var(cdt(x),cdt(y),...))
}
var.rcomp <- var.acomp
var.aplus <- var.acomp
var.rplus <- var.acomp

var.rmult <- function(x,y=NULL,...) {
  return(var(unclass(x),unclass(y),...))
  # XXXXXXXXX
  stop("Not yet implemented")
  if( is.null(y) ) {
    # Variance
    


  } else {
    # Covariance

  }

}

cov         <- function(x,y=x,...) UseMethod("cov",x)
cov.default <- stats::cov
formals(cov.default) <- c(formals(cov.default),alist(...= ))
cov.acomp   <- var.acomp
cov.rcomp <- cov.acomp
cov.aplus <- cov.acomp
cov.rplus <- cov.acomp
cov.rmult <- var.rmult

cor <- function(x,y=NULL,...) UseMethod("cor",x)
cor.default <- stats::cor
formals(cor.default) <- c(formals(cor.default),alist(...= ))

cor.acomp <- function(x,y=NULL,...) {
  if( has.missings(x) || has.missings(y) )
    cor(cdt(x),cdt(y),...)
  else
    cor(cdt(x),cdt(y),...)
}

cor.rcomp <- cor.acomp 
cor.aplus <- cor.acomp
cor.rplus <- cor.acomp
cor.rmult <- function(x,y=NULL,...) {
  cor(unclass(x),unclass(cdt(y)),...)
  
}


powerofpsdmatrix <- function(M,p,...) {
  s <- svd(M,...)
  d <- ifelse( abs(s$d)>max(abs(s$d))*1E-10, s$d^p,0)
  s$u %*% gsi.diagGenerate(d) %*% t(s$v)
}

mvar <- function(x,...) UseMethod("mvar",x)
mcov <- function(x,...) UseMethod("mcov",x)
mcor <- function(x,...) UseMethod("mcor",x)
msd  <- function(x,...) UseMethod("msd",x)

mvar.default <- function(x,y=NULL,...) {
  sum(gsi.diagExtract(var(x,y,...)))
}

mcov.default <- function(x,y=x,...) {
  sum(abs(svd(cov(idt(x),idt(y),...))$d))
}

msd.default <- function(x,y=NULL,...) {
  sqrt(mean(gsi.diagExtract(var(idt(x),y=NULL,...))))
}

mcor.default <- function(x,y,...) {
  ix <- scale(idt(x),center=TRUE,scale=FALSE)
  ix <- ix %*% powerofpsdmatrix(var(ix),-1/2)
  iy <- scale(idt(y),center=TRUE,scale=FALSE)
  iy <- iy %*% powerofpsdmatrix(var(iy),-1/2)
  mcov(ix,iy)
}



summary.acomp <- function( object,... ) {
  W <- clo(gsi.plain(object))
  Wq <- apply(gsi.recodeM2C(W,BDL=NaN,SZ=NaN,MAR=NaN,MNAR=NaN),
              1,function(w) outer(w,w,"/"))
  dim(Wq)<-c(ncol(W),ncol(W),nrow(W))
  dimnames(Wq) <- list(colnames(W),colnames(W),NULL)
  
  structure(list(mean=mean(acomp(W)),
       mean.ratio=apply(Wq,1:2,function(x) exp(mean(log(x)))),
       variation=variation.acomp(acomp(W)),
       expsd=exp(sqrt(variation.acomp(acomp(W)))),
       min=apply(Wq,1:2,min),
       q1 =apply(Wq,1:2,quantile,probs=0.25),
       med=apply(Wq,1:2,median),
       q3 =apply(Wq,1:2,quantile,probs=0.75),
       max=apply(Wq,1:2,max)
       ),class="summary.acomp")
       
}

summary.aplus <- function( object,...,digits=max(3, getOption("digits")-3)  ) {
  object <- ilt(object)
  erg <- sapply(data.frame(object),summary,...,digits=18)
  erg <- apply(erg,1:2,exp)
  erg <- apply(erg,1:2,signif,digits=digits)
  class(erg) <- c("summary.aplus",class(erg))
  erg       
}

summary.rplus <- function( object,...  ) {
  object <- iit(object)
  erg <- sapply(data.frame(object),summary,...)
  class(erg) <- c("summary.rplus",class(erg))
  erg       
}

summary.rmult <- function( object,...  ) {
  object <- unclass(object)
  erg <- sapply(data.frame(object),summary,...)
  class(erg) <- c("summary.rmult",class(erg))
  erg       
}

summary.rcomp <- function( object,...) {
  object <- clo(gsi.plain(object)) 
  erg <- sapply(data.frame(object),summary,...)
  class(erg) <- c("summary.rcomp",class(erg))
  erg       
}



vp.logboxplot <- function(x,y,...,dots=FALSE,boxes=TRUE,xlim,ylim,log=TRUE,notch=FALSE) {
    if( boxes ) {
      stats <- boxplot(split(log(y),x),plot=FALSE)
      stats$stats <- exp(stats$stats)
      stats$conf  <- exp(stats$conf)
      stats$out   <- exp(stats$out)
      bxp(stats,add=TRUE,notch=notch)
    }
    if( dots  ) points(x,y,...)
}



vp.boxplot <- function(x,y,...,dots=FALSE,boxes=TRUE,xlim,ylim,log,notch=FALSE) {
    if( boxes ) boxplot(split(y,x),add=TRUE,notch=notch)
    if( dots  ) points(x,y,...)
}


gsi.textpanel <- function(x,y,lab,...) {
  par(usr=c(0,1,0,1),xlog=FALSE,ylog=FALSE)
  text(0.5,0.5,lab,...)
}

boxplot.acomp <- function(x,fak=NULL,...,
                          xlim=x.lim,ylim=c(minq,maxq),log=TRUE,
                          panel=vp.logboxplot,dots=!boxes,boxes=TRUE) {
  X <- acomp(x)
  if( is.null(fak) )
    fak <- factor(rep("",nrow(X)))
  if( is.factor(fak) )
    x.lim <- c(0,nlevels(fak)+1)
  else {
    x.lim <- range(fak)
    boxes <- F
    dots  <- T
  }
  if( is.function(panel) )
    panel <- list(panel)
  ipanel <- function(x,y,...) {
    a <- gsi.recodeM2Clean(unclass(X)[,gsi.mapfrom01(log(x))])
    b <- gsi.recodeM2Clean(unclass(X)[,gsi.mapfrom01(log(y))])
    for( thispanel in panel ) 
      thispanel(fak,b/a,...,dots=dots,boxes=boxes)
  }
  su <- summary.acomp(X)
  minq <- min(su$min)
  maxq <- max(su$max)
  mm <- exp(sapply(1:NCOL(X),gsi.mapin01))
  colnames(mm) <- colnames(X)
  ipairs <- function (x, labels, panel = ipanel, ..., 
                      font.main = par("font.main"),
                      cex.main = par("cex.main"), diag.panel = NULL, 
                      text.panel = textPanel,
                      label.pos = 0.5 , cex.labels = NULL, 
                      font.labels = 1, gap = 1,xlim,ylim,log) {
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                     y, txt, cex = cex, font = font)
    localAxis <- function(side, xpd, bg, ...) axis(side, xpd = NA, 
                                                   ...)
    panel <- match.fun(panel)
    nc <- ncol(x)
    labels <- colnames(x)
    if (is.null(labels)) 
      labels <- paste("var", 1:nc)
    has.labs <- ! is.null(labels)
    oma <- c(4, 4, 4, 4)
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    for (i in 1:nc ) for (j in 1:nc) {
      plot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
           type = "n", ...,xlim=xlim,ylim=ylim,log=log)
      box()
      mfg <- par("mfg")
      if (i == j) {
        if (has.labs) {
          par(usr = c(0, 1, 0, 1))
          if (is.null(cex.labels)) {
            l.wid <- strwidth(labels, "user")
            cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
          }
          text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                     font = font.labels)
        }
      }
      else  
        panel(as.vector(x[, j]), as.vector(x[, i]), ...)
      if (any(par("mfg") != mfg)) 
        stop("The panel function made a new plot")
    }
    invisible(NULL)
  }
  
  ipairs(mm,labels=colnames(X),panel=ipanel,...,log=ifelse(log,"y",""),ylim=ylim,xlim=xlim,text.panel=gsi.textpanel)
  
  
}

boxplot.rcomp <- function(x,fak=NULL,...,
                         xlim=x.lim,ylim=c(0,1),log=FALSE,panel=vp.boxplot,dots=!boxes,boxes=TRUE) {
  X <- acomp(x)
  if( is.null(fak) )
    fak <- factor(rep("",nrow(X)))
  if( is.factor(fak) )
    x.lim <- c(0,nlevels(fak)+1)
  else {
    x.lim <- range(fak)
    boxes <- F
    dots  <- T
  }
  if( is.function(panel) )
    panel <- list(panel)
  if( missing(ylim) && log )
    ylim <- c(minq,maxq)
  ipanel <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(log(x))]
    b <- unclass(X)[,gsi.mapfrom01(log(y))]
    for( thispanel in panel ) 
      thispanel(fak,if(log) log(b/(a+b)) else b/(a+b),
                ...,dots=dots,boxes=boxes)
  }
  su <- summary.acomp(X)
  minq <- log(min(su$min)/(min(su$min)+max(su$max)))
  maxq <- log(max(su$max)/(min(su$min)+max(su$max)))
  
  ipairs <- function (x, labels, panel = points, ..., 
                      font.main = par("font.main"),
                      cex.main = par("cex.main"), diag.panel = NULL, 
                      text.panel = textPanel,
                      label.pos = 0.5 , cex.labels = NULL, 
                      font.labels = 1, gap = 1,xlim,ylim,log="") {
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                     y, txt, cex = cex, font = font)
    localAxis <- function(side, xpd, bg, ...) axis(side, xpd = NA, 
                                                   ...)
    panel <- match.fun(panel)
    nc <- ncol(x)
    labels <- colnames(x)
    if (is.null(labels)) 
      labels <- paste("var", 1:nc)
    has.labs <- ! is.null(labels)
    oma <- c(4, 4, 4, 4)
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    for (i in 1:nc ) for (j in 1:nc) {
      plot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
           type = "n", ...,xlim=xlim,ylim=ylim,log=log)
      box()
      mfg <- par("mfg")
      if (i == j) {
        if (has.labs) {
          par(usr = c(0, 1, 0, 1))
          if (is.null(cex.labels)) {
            l.wid <- strwidth(labels, "user")
            cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
          }
          text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                     font = font.labels)
        }
      }
      else  
        panel(as.vector(x[, j]), as.vector(x[, i]), ...)
      if (any(par("mfg") != mfg)) 
        stop("The panel function made a new plot")
    }
    invisible(NULL)
  }


  ipairs(exp(sapply(1:NCOL(X),gsi.mapin01)),labels=colnames(X),panel=ipanel,...,ylim=ylim,xlim=xlim)
  

}

boxplot.rplus <- function(x,fak=NULL,...,ylim=c(0,max(x)),log=FALSE) {
  if( !is.null(fak) )
    warning("Spliting not yet implemente in boxplot.rplus")
  boxplot(as.data.frame(x),...,ylim=ylim,log=if(log) "y" else "")
}

boxplot.aplus <- function(x,fak=NULL,...,log=TRUE) {
  if( !is.null(fak) )
    warning("Spliting not yet implemente in boxplot.aplus")
  stats <- boxplot(as.data.frame(ilt(x)),plot=FALSE)
  delog <- function(x) {if(is.list(x)) lapply(x,delog) else exp(x)}
  stats$stats <- exp(stats$stats)
  stats$conf <- exp(stats$conf)
  stats$out <- exp(stats$out)
  invisible(bxp(stats,...,log=if(log) "y" else ""))
}



vp.qqnorm <- function(x,y,...,alpha=NULL) {
  usr <- par("usr")
  usr[1:2] <- range(qnorm(ppoints(length(y))))
  usr[3:4] <- range(y)
  par( usr=usr )
  if( !is.null(alpha) && is.factor(x) ) 
    alpha <- alpha/nlevels(x)
  reject <- FALSE
  if( is.factor(x)) {
    for( k in split(y,x) ) {
      if( !is.null(alpha) && shapiro.test(k)$p < alpha )
        reject<-TRUE
      lines(qnorm(ppoints(length(k))),sort(k),...)
    }
  } else { 
    if( !is.null(alpha) && shapiro.test(y)$p < alpha )
        reject<-TRUE
    points(qnorm(ppoints(length(y))),sort(y),...)
  }
  qqline(y)
  if( reject )
    title(main="!",col.main="red")
    
}

qqnorm.acomp <- function(y,fak=NULL,...,panel=vp.qqnorm,alpha=NULL) {
  X <- acomp(y)
  if( !is.null(alpha) )
    alpha <- alpha/((nrow(X)*(nrow(X)-1)/2))
  if( is.function(panel) )
    panel <- list(panel)
  ipanel <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    v <- log(b/a)
    for( thispanel in panel )
      thispanel(fak,v,...,alpha=alpha)
  }
  pairs(sapply(1:NCOL(X),gsi.mapin01),labels=colnames(X),panel=ipanel,...)
}

qqnorm.aplus <- function(y,fak=NULL,...,panel=vp.qqnorm,alpha=NULL) {
  X <- aplus(y)
  if( is.function(panel) )
    panel <- list(panel)
  if( !is.null(alpha) )
    alpha <- alpha/(nrow(X)^2)
  ipanelupper <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    for( thispanel in panel )
      thispanel(fak,log(b/a),...,alpha=alpha)
  }
  ipanellower <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    for( thispanel in panel )
      thispanel(fak,log(b*a),...,alpha=alpha)
  }
  ipaneldiag <- function(x,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    for( thispanel in panel )
      thispanel(fak,log(a),...,alpha=alpha)
  }
  itextpanel <- function(x,y,lab,...) {
    par(usr=c(0,1,0,1),xlog=FALSE,ylog=FALSE)
    text(0.1,0.9,lab,adj=c(0,1),...)
  }

  pairs(sapply(1:NCOL(X),gsi.mapin01),labels=colnames(X),lower.panel=ipanellower,upper.panel=ipanelupper,diag.panel=ipaneldiag,text.panel=itextpanel,...)
}



qqnorm.rcomp <- function(y,fak=NULL,...,panel=vp.qqnorm,alpha=NULL) {
  X <- rcomp(y)
  if( is.function(panel) )
    panel <- list(panel)
  if( !is.null(alpha) )
    alpha <- alpha/(nrow(X)^2)
  ipanelupper <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    for( thispanel in panel )
      thispanel(fak,b-a,...,alpha=alpha)
  }
  ipanellower <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    for( thispanel in panel )
      thispanel(fak,b+a,...,alpha=alpha)
  }
  ipaneldiag <- function(x,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    for( thispanel in panel )
      thispanel(fak,a,...,alpha=alpha)
  }
  itextpanel <- function(x,y,lab,...) {
    par(usr=c(0,1,0,1),xlog=FALSE,ylog=FALSE)
    text(0.1,0.9,lab,adj=c(0,1),...)
  }

  pairs(sapply(1:NCOL(X),gsi.mapin01),labels=colnames(X),lower.panel=ipanellower,upper.panel=ipanelupper,diag.panel=ipaneldiag,text.panel=itextpanel,...)
}

qqnorm.rplus <- function(y,fak=NULL,...,panel=vp.qqnorm,alpha=NULL) {
  X <- rplus(y)
  if( is.function(panel) )
    panel <- list(panel)
  if( !is.null(alpha) )
    alpha <- alpha/(nrow(X)^2)
  ipanelupper <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    for( thispanel in panel )
      thispanel(fak,b-a,...,alpha=alpha)
  }
  ipanellower <- function(x,y,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    b <- unclass(X)[,gsi.mapfrom01(y)]
    for( thispanel in panel )
      thispanel(fak,b+a,...,alpha=alpha)
  }
  ipaneldiag <- function(x,...) {
    a <- unclass(X)[,gsi.mapfrom01(x)]
    for( thispanel in panel )
      thispanel(fak,a,...,alpha=alpha)
  }
  itextpanel <- function(x,y,lab,...) {
    par(usr=c(0,1,0,1),xlog=FALSE,ylog=FALSE)
    text(0.1,0.9,lab,adj=c(0,1),...)
  }

  pairs(sapply(1:NCOL(X),gsi.mapin01),labels=colnames(X),lower.panel=ipanellower,upper.panel=ipanelupper,diag.panel=ipaneldiag,text.panel=itextpanel,...)
}



gsi.drop  <-function(X,drop) if( drop ) drop(X) else X

is.acomp <- function(x) inherits(x,"acomp")

is.rcomp <- function(x) inherits(x,"rcomp")

is.aplus <- function(x) inherits(x,"aplus")

is.rplus <- function(x) inherits(x,"rplus")
   
is.rmult <- function(x) inherits(x,"rmult")


perturbe <- function( x,y ) {
  acomp(gsi.mul(x,y))
}

perturbe.aplus <- function(x,y) {
  aplus(gsi.mul(x,y))
}



gsi.add <- function( x,y ) {
  if( length(dim(x)) == 2 )
    if( length(dim(y)) == 2 )
      unclass(x)+unclass(y)
    else
      unclass(x)+rep(c(y),rep(NROW(x),length(y)))
  else if( length(dim(y)) == 2 )
      unclass(y)+rep(c(x),rep(NROW(y),length(x)))
  else
    unclass(x)+unclass(y)
}

gsi.sub <- function( x,y ) {
 # drop <- length(dim(x)) < 2 && length(dim(y)) < 2
  if( length(dim(x)) == 2 )
    if( length(dim(y)) == 2 )
      unclass(x)-unclass(y)
    else
      unclass(x)-rep(c(y),rep(NROW(x),length(y)))
  else if( length(dim(y)) == 2 )
      unclass(y)-rep(c(x),rep(NROW(y),length(x)))
  else
    unclass(x)-unclass(y)
}

gsi.mul <- function( x,y ) {
  if( length(dim(x)) == 2 )
    if( length(dim(y)) == 2 )
      unclass(x)*unclass(y)
    else
      unclass(x)*rep(c(y),rep(NROW(x),length(y)))
  else if( length(dim(y)) == 2 )
      unclass(y)*rep(c(x),rep(NROW(y),length(x)))
  else
    unclass(x)*unclass(y)
}

gsi.div <- function( x,y ) {
  if( length(dim(x)) == 2 )
    if( length(dim(y)) == 2 )
      unclass(x)/unclass(y)
    else
      unclass(x)/rep(c(y),rep(NROW(x),length(y)))
  else if( length(dim(y)) == 2 )
      unclass(y)/rep(c(x),rep(NROW(y),length(x)))
  else
    unclass(x)/unclass(y)
}


power.acomp <- function(x,s) {
  if( is.acomp(s) || is.rcomp(s))
    stop("power.acomp is scalar multiplication only")
  if( !is.matrix(x) || nrow(x) ==1 ) {
    if( length(s)>1 )
      x <- matrix(x,byrow=TRUE,ncol=length(x),nrow=length(s))
  } else {
    if( length(s) > 1 && length(s)!= nrow(x) )
      warning("lengths don't match in power.acomp")
  }
  acomp(unclass(x)^c(s)) 
}


"+.acomp" <- function(x,y) {
  acomp(gsi.mul(x,y))
}

"-.acomp" <- function(x,y) {
  if( missing(y) )
    acomp(1/unclass(x))
  else 
    acomp(gsi.div(x,y))
}

"*.acomp" <- function(x,y) {
  if( is.acomp(x) && !is.acomp(y) )
    power.acomp(x,y)
  else if( is.acomp(y)&& !is.acomp(x) )
    power.acomp(y,x)
  else
    stop("the powertransform performed in *.acomp only operates on acomps and scalar")
}

"/.acomp" <- function(x,y) {
  if( is.acomp(x) && !is.acomp(y) )
    power.acomp(x,1/unclass(y))
  else
    stop("/.acomp only operates on acomp / numeric")
}

"+.aplus" <- function(x,y) {
    aplus(gsi.mul(x,y))
}

"-.aplus" <- function(x,y) {
  if( missing(y) )
    return(aplus(1/unclass(y)))
  else
    aplus( gsi.div(x,y) )
}

"*.aplus" <- function(x,y) {
  if( is.aplus(x)&& !is.aplus(y) )
    power.aplus(x,y)
  else if( is.aplus(y)&& !is.aplus(x) )
    power.aplus(y,x)
  else
    stop("*.aplus only operates on aplus and scalar")
}

"/.aplus" <- function(x,y) {
  if( is.aplus(x) && !is.aplus(y) )
    power.aplus(x,1/unclass(y))
  else
    stop("/.aplus only operates on aplus and scalar")
}


"+.rcomp" <- function(x,y) {
#  warning("+ is meaningless for rcomp")
  if( is.rcomp(x) )
    if( is.rcomp(y) ) {
      stop("+ is meaningless for two rcomp objects")
   } else if( is.rcomp(x)) {
      rmult(clo(x))+rmult(y)
    } else if( is.rcomp(y) ) {
      rmult(x)+rmult(clo(y))
    } else
      stop("Why are we here in +.rcomp without rcomp?")
  #rcomp(gsi.add(x,y))
}

"-.rcomp" <- function(x,y) {
  if( missing(y) )
    rmult(-unclass(x))
  else
    rmult(gsi.sub(x,y))
}

"*.rcomp" <- function(x,y) {
  if( is.rcomp(x) && is.rcomp(y) )
    rcomp(gsi.mul(x,y))
  else if( is.rcomp(x) )
    rplus(x)*y
  else if( is.rcomp(y) )
    rplus(y)*x
  else
    stop("undefined combination of arguments for *.rcomp")
}

"/.rcomp" <- function(x,y) {
  if( is.rcomp(x) && is.rcomp(y) )
    rcomp(gsi.div(x,y))
  else if( is.rcomp(x) )
    rplus(x)/y
  else
    stop("undefined combination of arguments for /.rcomp")
}

"+.rplus" <- function(x,y) {
  if( is.rplus(x) && is.rplus(y) )
    rplus(gsi.add(x,y))
  else
    rmult(gsi.add(x,y))
}

"-.rplus" <- function(x,y) {
  if( missing(y) )
    rmult(-unclass(x))
  else
    rmult(gsi.sub(x,y))
}


"*.rplus" <- function(x,y) {
  if( is.rplus(x) && is.rplus(y) )
    rplus(gsi.mul(x,y))
  else if( is.rplus(x) )
    mul.rplus(x,y)
  else if( is.rplus(y) )
    mul.rplus(y,x)
  else
    stop("undefined combination of arguments for *.rplus")
}

"/.rplus" <- function(x,y) {
  if( is.rplus(x) && is.rplus(y) )
    rplus(gsi.div(x,y))
  else if( is.rcomp(x) )
    mul.rplus(rplus(x),1/unclass(y))
  else
    stop("undefined combination of arguments for /.rcomp")
}

"+.rmult" <- function(x,y) {
  rmult(gsi.add(x,y))
}

"-.rmult" <- function(x,y) {
  if( missing(y) )
    rmult(-unclass(x))
  else
    rmult(gsi.sub(x,y))
}


"*.rmult" <- function(x,y) {
  if( is.rmult(x) && is.rmult(y) )
    rmult(gsi.mul(x,y))
  else
    rmult(unclass(x)*unclass(y))
}

"/.rmult" <- function(x,y) {
  if( is.rmult(x) && is.rmult(y) )
    rmult(gsi.div(x,y))
  else 
    rmult(unclass(x)/unclass(y))
}

"%*%" <- function(x,y) UseMethod("%*%",structure(c(),class=c(class(x),class(y))))


#gsi.internaltmp <- get("%*%",pos="package:base")
#formals(gsi.internaltmp) <- formals(get("%*%"))
#"%*%.default" <- gsi.internaltmp

"%*%.default" <- function(x,y) base::"%*%"(x,y)

"%*%.rmult" <- function(x,y) {
  if( is.rmult(y) )
    if( is.rmult(x) ) 
      c(gsi.mul(x,y) %*% rep(1,gsi.getD(x)))
    else if( is.matrix(x) ) 
      rmult(gsi.simshape(oneOrDataset(y) %*% t(x),y))
    else
      c(oneOrDataset(y) %*% x) 
  else if( is.matrix(y) )
      rmult(gsi.simshape(oneOrDataset(x) %*% y,x))
  else
      c(oneOrDataset(x) %*% y) 
  }

"%*%.acomp" <- function(x,y) {
  if( is.acomp(y) )
    if( is.acomp(x) ) 
      cdt(x) %*% cdt(y)
    else if( is.matrix(x) ) {
      if( nrow(x) == gsi.getD(y) )
        clr.inv(x %*% clr(y))
      else
        ilr.inv(x %*% ilr(y))
    }
    else
      stop( "%*%.acomp is only defined for special combinations I" )
  else if( is.acomp(x) ) {
    if( is.matrix(y) ) {
      if( ncol(y) == gsi.getD(x) )
        clr.inv(clr(x) %*% y )
      else
        ilr.inv(ilr(x) %*% y )
    }
  else
      stop( "%*%.acomp is only defined for special combinations II" )
  }
  else
      stop( "%*%.acomp is only defined for special combinations III" )
    
}

"%*%.aplus" <- function(x,y) {
  if( is.aplus(y) )
    if( is.aplus(x) ) 
      cdt(x) %*% cdt(y)
    else if( is.matrix(x) ) {
        ilt.inv(x %*% ilt(y))
    }
    else
      stop( "%*%.acomp is only defined for special combinations I" )
  else if( is.aplus(x) ) {
    if( is.matrix(y) ) {
        ilt.inv(ilt(x) %*% y )
    }
  else
      stop( "%*%.aplus is only defined for special combinations II" )
  }
  else
      stop( "%*%.aplus is only defined for special combinations III" )
    
}


convex.rcomp <- function(x,y,alpha=0.5) {
  rcomp( alpha*x + (1-alpha)*y )
}


mul.rplus <- function(x,r) {
  if( all(r>=0) )
    rplus(unclass(x)*r)
  else
    rmult(unclass(x)*r)
}

power.aplus <- function(x,r) {
  aplus(unclass(x)^r) 
}


gsi.expandrcomp <- function(x,alpha) {
  cpt.inv(cpt(x)*alpha)
}

endmemberCoordinates <- function(X,...) UseMethod("endmemberCoordinates")

endmemberCoordinates.default <- function(X,endmembers=diag(gsi.getD(X)),...) {
  class(endmembers) <- class(X)
  X <- oneOrDataset(idt(X))
  A <- t(unclass(idt(endmembers)))
  erg <- solve( rbind(cbind(t(A)%*%A,1),c(rep(1,ncol(A)),0)),
               rbind(t(A)%*%t(unclass(X)),1))
  erg <- rmult(t(erg[-nrow(erg),,drop=FALSE]))
  colnames(erg) <- rownames(endmembers)
  erg
}

endmemberCoordinates.acomp <- function(X,endmembers=clr.inv(diag(gsi.getD(X))),...) {
  ep <- ilr(endmembers)
  rownames(ep) <- rownames(endmembers)
  endmemberCoordinates(idt(X),ep,...)
}

endmemberCoordinates.aplus <- function(X,endmembers,...) {
  ep <- ilt(endmembers)
  rownames(ep) <- rownames(endmembers)
  endmemberCoordinates(idt(X),ep,...)
}


endmemberCoordinates.rplus <- function(X,endmembers,...) {
  ep <- iit(endmembers)
  rownames(ep) <- rownames(endmembers)
  endmemberCoordinates(idt(X),ep,...)
}


endmemberCoordinatesInv <- function(K,endmembers,...) UseMethod("endmemberCoordinatesInv",endmembers)

endmemberCoordinatesInv.rmult <- function(K,endmembers,...) {
  rmult(t(t(unclass(endmembers)) %*% t(unclass(K))))
}

endmemberCoordinatesInv.acomp <- function(K,endmembers,...) {
  ilr.inv(endmemberCoordinatesInv(K,ilr(endmembers)))
}

endmemberCoordinatesInv.rcomp <- function(K,endmembers,...) {
  ipt.inv(endmemberCoordinatesInv(K,ipt(endmembers)))
}


endmemberCoordinatesInv.aplus <- function(K,endmembers,...) {
  ilt.inv(endmemberCoordinatesInv(K,ilt(endmembers)))
}

endmemberCoordinatesInv.rplus <- function(K,endmembers,...) {
  iit.inv(endmemberCoordinatesInv(K,iit(endmembers)))
}



scale.acomp <- function( x,center=TRUE, scale=TRUE ) {
  W <- x
  if( center ) {
    W <- clr.inv( scale(clr(W),center=center,scale=FALSE) )
    if( scale )
      W <- power.acomp(W,as.numeric(scale)/
                       sqrt(sum(gsi.diagExtract(var(clr(W))))))
  } else if( scale ) {
    mean <- c(mean(acomp(W)))
    W <- perturbe(power.acomp(perturbe(W,1/mean),as.numeric(scale)/sqrt(sum(gsi.diagExtract(var(clr(W)))))),mean)
  }
  W
}

scale.rcomp <- function( x,center=TRUE, scale=TRUE ) {
  W <- x
  if( center ) {
    W <- cpt.inv( scale(cpt(W),center=center,scale=FALSE) )
    if( scale )
      W <- gsi.expandrcomp(W,as.numeric(scale)/sqrt(sum(gsi.diagExtract(var(cpt(W))))))
  } else if( scale ) {
    mean <- c(mean(rcomp(W)))
    W <- gsi.add(mean,gsi.sub(W,mean)/sqrt(sum(gsi.diagExtract(var(cpt(W))))))
  }
  W
}

scale.aplus <- function( x,center=TRUE, scale=TRUE ) {
  W <- x
  if( center ) {
    W <- ilt.inv( scale(ilt(W),center=center,scale=FALSE) )
    if( scale )
      W <- power.aplus(W,as.numeric(scale)/sqrt(sum(gsi.diagExtract(var(ilt(W))))))
  } else if( scale ) {
    mean <- c(mean(aplus(W)))
    W <- perturbe.aplus(power.aplus(perturbe.aplus(W,1/mean),as.numeric(scale)/sqrt(sum(gsi.diagExtract(var(ilt(W)))))),mean)
  }
  W
}

scale.rplus <- function( x,center=TRUE, scale=TRUE ) {
   rmult(scale(gsi.plain(x),center=center,scale=scale))
}

scale.rmult <- function( x,center=TRUE, scale=TRUE ) {
   rmult(scale(gsi.plain(x),center=center,scale=scale))
}

normalize <- function(x,...) UseMethod("normalize",x)
normalize.default <- function(x,...) x/norm(x)

if( !exists("norm")) norm <- function(x,...) UseMethod("norm",x)

norm.default <- function(X,...) {
  sqrt( sum(X^2) )
}

norm.acomp <- function(X,...) {
  norm.rmult(cdt(X),...)
}
norm.rcomp <- norm.acomp
norm.aplus <- norm.acomp
norm.rplus <- norm.acomp
norm.rmult <- function(X,...) {
   sqrt(X %*% X)
}

dist <- function(x,...) UseMethod("dist")
dist.default <- function(x,...) stats::dist(cdt(x),...)


scalar <- function(x,y) UseMethod("scalar")

scalar.default <- function(x,y) {
  x <- cdt(x)
  y <- cdt(y)
  tmp <- gsi.mul(oneOrDataset(x,y), oneOrDataset(y,x)) 
  c( tmp %*% rep(1,NCOL(tmp)))
}


clr <- function( x ,... ) {
  W <- oneOrDataset(x)
  nmv <- is.NMV(W)
  LOG <- unclass(log(ifelse(nmv,W,1)))
  erg <- ifelse(nmv,LOG-gsi.rsum(LOG)/gsi.rsum(nmv),0)
  #M   <- gsi.rsum(gsi.recodeM2C(W,LOG,BDL=0,SZ=0,MAR=0,MNAR=0))
  rmult(gsi.simshape(erg,x),orig=x) 
}

clr.inv <- function( z,... ) {
  acomp( gsi.recodeC2M(exp(z),ninf=BDLvalue,nan=MARvalue,na=MNARvalue) )
}

ult <- function( x,... ) {
  ilt(clo(x),...)
}

ult.inv <- clr.inv

Kappa <- function( x, ...) {
  W <- oneOrDataset(x)
  (clr(W)-ult(W))[,1]
}

#gsi.ilrBase <- function(D) {
#  if( D==1 )
#    return(matrix(nrow=0,ncol=0))
#  tmp <- diag(D) - 1/(D)* matrix(1,ncol=D,nrow=D)
#  for( i in 1:(NCOL(tmp)-1)  ) {
#    tmp[,i] <- tmp[,i]/sqrt(sum(tmp[,i]^2))
#    rest <- (i+1):NCOL(tmp)
#    if( length(rest) != 1 ) {
#      tmp[,rest]  <-tmp[,rest,drop=FALSE] - tmp[,rep(i,length(rest)),drop=FALSE]%*%
#        gsi.diagGenerate( c(t(tmp[,i])%*%tmp[,rest,drop=FALSE] ) )
#    } 
#  }
#  tmp[,-NROW(tmp)]
#}

gsi.ilrBase <-function(D){
  if(D<=1)
    return(matrix(nrow=0,ncol=0))
  else 
    t(unclass(normalize(rmult(t(contr.helmert(n=D))))))
}




                                        # Old ilrBase
#ilrBase <- function( x=NULL , z=NULL , D = NULL ) {
#  if( missing(D) )
#    D <- if(is.null(x))
#      NCOL(oneOrDataset(z))+1
#    else
#      NCOL(oneOrDataset(x))
#  while( D > length(ilrBaseList) )
#    ilrBaseList <<- c(ilrBaseList,gsi.ilrBase(length(ilrBaseList)+1))
#  ilrBaseList[[D]]gsi.OrderIlr
#}
# ilrBaseList <- lapply(1:20,gsi.ilrBase)


ilrBase <- function(x=NULL, z=NULL, D=NULL, method="basic"){
 if (method=="basic"){
    if (missing(D))
        D <- if (is.null(x))
            NCOL(oneOrDataset(z)) + 1
        else NCOL(oneOrDataset(x))
#     while (D > length(ilrBaseList)) ilrBaseList <<- c(ilrBaseList,
#         gsi.ilrBase(length(ilrBaseList) + 1))
#     V=ilrBaseList[[D]]
  V = gsi.ilrBase(D)
 } #end if basic
 if (method=="balanced"){
    # build the merge structure (as in hclust)
    if(is.null(D)==TRUE){D=max(ncol(x),ncol(z))}
    M = c(-c(1:D),c(floor(D/2):1),c((floor(D/2)+1):(D-2)))
    M = matrix(M, byrow=TRUE, nrow=D-1, ncol=2)
    V = gsi.merge2signary(M)
    V = gsi.buildilrBase(V)
 }#end if balanced
 if (method=="optimal"){
  if(length(dim(x))<2){
   warning("method 'optimal' requires a data set")
  }
  else{
    V = gsi.optimalilrBase(x)
    V = gsi.buildilrBase(V)

  }
 } #end if optimal
 return(V)
}

gsi.buildilrBase <-function(W=c(1,-1)){
  # builds an ilr base from a matrix of 1, -1 and 0 (a partition)
  if(length(W)<2){
    return(ilrBase(D=1))
  }
if(length(dim(W))==0){
  return(ilrBase(D=2))
  }
  if(length(dim(W))>0){
   W = as.matrix(W)
   nc = ncol(W)
   D = nrow(W)
   isPos = (W>0)
   isNeg = (W<0)
   nPos = matrix(1,D,D) %*% isPos
   nNeg = matrix(1,D,D) %*% isNeg
   W = (isPos * nNeg - isNeg * nPos)
   nn = sapply(1:nc,function(i){1/norm.rmult(W[,i])})
   nn = matrix(nn, ncol=ncol(W), nrow=nrow(W), byrow=TRUE)
   W = W * nn
   return(W)
  }
}
gsi.signary2ilrBase <- gsi.buildilrBase


gsi.optimalilrBase <- function(x){
  # Takes as data a binary data frame: 0=case unobserved, 1=case observed
  # performs a cluster analysis of variables on it
  # and recodes the structure to a signary matrix (pre-ilr)
  if(is.null(attr(x,"Losts"))==TRUE){
    Ones = !oneOrDataset(is.infinite(log(as.matrix(x))),x) # find zeroes
  }
  if(is.null(attr(x,"Losts"))==FALSE){
    Ones = !attr(x,"Losts")
  }
  h= hclust(dist(t(Ones)))
  V = gsi.merge2signary(h$merge)
  return(V)
}


gsi.merge2signary <- function(M){
  # takes a merge structure (as explained in hclust) and converts it to a sign matrix (encoding the ilr partition: 0=no influence, 1=numerator, -1=denominator)
 V = matrix(0,ncol=nrow(M)+1,nrow=nrow(M))
  for(i in 1:nrow(M)){
   for(j in 1:2){
    weight = (-1)^j
    k = M[i,j]
    if(k<0){
      V[i,abs(k)] = weight
    } # for singletons
    if(k>0){ # for groups
      take = as.logical(V[k,])
      V[i,take] = rep(weight,sum(take))
    }
  }
 }
return(t(V))
}




ilr    <- function( x , V=ilrBase(x),... ) {
  rmult(clr2ilr( clr(oneOrDataset(x)),V ))
}

ilr.inv <- function( z, V=ilrBase(z=z),...,orig=NULL) {
  erg <- clr.inv( ilr2clr(z,V) )
  if( ! is.null(orig) && gsi.getD(erg) == gsi.getD(orig) ) {
    names(erg)<-names(orig)
  }
  erg
}

alr <- function( x ,...) {
W <- unclass(clo(oneOrDataset(x)))
V <- gsi.recodeM2C(W,log(W),BDL=-Inf,SZ=NaN,MAR=NaN,MNAR=NA)
rmult(gsi.simshape( V[,-NCOL(V),drop=FALSE] - c(V[,NCOL(V)]), x))
}

alr.inv <- function( z ,...,orig=NULL) {
  Z <- cbind(oneOrDataset(z),0)
  erg <- acomp(gsi.simshape( clo(gsi.recodeC2M(Z,exp(Z),
                                        ninf=BDLvalue,
                                        inf =MNARvalue,
                                        nan =MARvalue,
                                        na  =MNARvalue
                                        )) , z ))
  if( ! is.null(orig) && gsi.getD(erg) == gsi.getD(orig) ) {
    names(erg)<-names(orig)
  }
  erg
}


apt <- function( x ,...) {
  W <- oneOrDataset(x)
  V <- gsi.recodeM2C(W,gsi.plain(clo( W )),BDL=0.0,SZ=0.0,MAR=NaN,MNAR=NA)
  rmult(gsi.simshape(V[,-NCOL(W)],x))
}

apt.inv <- function( z ,...,orig=NULL) {
  Z <- oneOrDataset(z)
  Z <- cbind(Z, 1 - gsi.recodeC2M(Z,na=0.0,nan=0.0) %*% rep(1,NCOL(Z)))
  erg <- rcomp(gsi.simshape( Z ,z ))
  if( ! is.null(orig) && gsi.getD(erg) == gsi.getD(orig) ) {
    names(erg)<-names(orig)
  }
  erg
}

cpt <- function( x ,...) {
  X <- oneOrDataset(x)
  rmult(clo(x) - 1/NCOL(X),orig=x)
}

cpt.inv <- function( z ,...) {
  if( abs(sum(z))>0.0001 )
    warning( "z not in cpt plane in cpt.inv")
  rcomp(z + 1/NCOL(oneOrDataset(z)))
}

ipt    <- function( x , V=ilrBase(x),...) {
  rmult(clr2ilr(cpt(x),V),orig=x)
  
}

ipt.inv <- function( z, V=ilrBase(z=z) ,...,orig=NULL) {
  erg<-cpt.inv( ilr2clr(z,V) )
  if( ! is.null(orig) && gsi.getD(erg) == gsi.getD(orig) ) {
    names(erg)<-names(orig)
  }
  erg
}

ucipt.inv <- function( z, V=ilrBase(z=z) ,...,orig=NULL) {
  tmp <- ilr2clr(z,V) + 1/(gsi.getD(z)+1)
  tmp[tmp<0]<-MNARvalue
  erg<-rcomp(tmp)
  if( ! is.null(orig) && gsi.getD(erg) == gsi.getD(orig) ) {
    names(erg)<-names(orig)
  }
  erg
}


ilt <- function( x ,...) {
  W <- gsi.plain(x)
  rmult(log(ifelse(is.NMV(W),W,1)),orig=x)
}

ilt.inv <- function( z ,...) {
  aplus(gsi.recodeC2M(z,exp(z),ninf=BDLvalue,nan=MARvalue,na=MNARvalue))
}

iit <- function( x ,...) {
  rmult( gsi.recodeM2C(x,BDL=0.0,SZ=0.0,MAR=0.0,MNAR=0.0 ) ,orig=x)
}

iit.inv <- function(z,...) {
  rplus(gsi.recodeC2M(z,na=MNARvalue,nan=MARvalue))
}

idt         <- function(x,...) UseMethod("idt",x)
idt.default <- function(x,...) x
idt.acomp   <- function(x,...) ilr(x,...) 
idt.rcomp   <- function(x,...) ipt(x,...) 
idt.aplus   <- ilt 
idt.rplus   <- iit 
idt.rmult   <- function(x,...) x
idt.factor  <- function(x,...) rmult(clr2ilr(cdt(x)))


cdt         <- function(x,...) UseMethod("cdt",x)
cdt.default <- function(x,...) x
cdt.acomp   <- clr 
cdt.rcomp   <- cpt 
cdt.aplus   <- ilt 
cdt.rplus   <- iit 
cdt.rmult   <- function(x,...) x
cdt.factor  <- function(x,...) {
  #x <- matrix(0,nrow=length(x),ncol=nlevels(x),dimnames=list(names(x),levels(x)))
  x[1:ncol(x)+unclass(x)] <- model.matrix(~-1+x)
  
  rmult(matrix(x,nrow=nrow(x),dimnames=dimnames(x)))
}

cdt.inv <- function(x,orig,...) UseMethod("cdt.inv",orig)
cdt.inv.default <- function(x,orig,...) x
cdt.inv.acomp   <- function(x,orig,...) clr.inv(x,...,orig=orig)
cdt.inv.rcomp   <- function(x,orig,...) cpt.inv(x,...,orig=orig)
cdt.inv.aplus   <- function(x,orig,...) ilt.inv(x,...,orig=orig)
cdt.inv.rplus   <- function(x,orig,...) iit.inv(x,...,orig=orig)
cdt.inv.rmult   <- function(x,orig,...) x

idt.inv <- function(x,orig,...) UseMethod("idt.inv",orig)
idt.inv.default <- function(x,orig,...) x
idt.inv.acomp   <- function(x,orig,...) ilr.inv(x,...,orig=orig)
idt.inv.rcomp   <- function(x,orig,...) ipt.inv(x,...,orig=orig)
idt.inv.aplus   <- function(x,orig,...) ilt.inv(x,...,orig=orig)
idt.inv.rplus   <- function(x,orig,...) iit.inv(x,...,orig=orig)
idt.inv.rmult   <- function(x,orig,...) x


variation <- function( x, ... ) UseMethod("variation",x)

variation.acomp <- function( x,... ) {
  co <-var(clr(x))
  d <- NCOL(x)
  va <-gsi.diagExtract(co)
  co1 <- matrix(rep(va,each=d),ncol=d)
  co2 <- matrix(rep(va,d),ncol=d)
  -2*co+co1+co2
  
}

variation.rcomp <- function( x ,...) {
  co <-var(cpt(x))
  d <- NCOL(x)
  va <-gsi.diagExtract(co)
  co1 <- matrix(rep(va,each=d),ncol=d)
  co2 <- matrix(rep(va,d),ncol=d)
  -2*co+co1+co2
  
}


variation.aplus <- function( x ,...) {
  co <-var(ilt(x))
  d <- NCOL(x)
  va <-gsi.diagExtract(co)
  co1 <- matrix(rep(va,each=d),ncol=d)
  co2 <- matrix(rep(va,d),ncol=d)
  -2*co+co1+co2
  
}

variation.rmult <- function( x ,...) {
  co <-var(iit(x))
  d <- NCOL(x)
  va <-gsi.diagExtract(co)
  co1 <- matrix(rep(va,each=d),ncol=d)
  co2 <- matrix(rep(va,d),ncol=d)
  -2*co+co1+co2
  
}
variation.rplus <- variation.rmult



is.BDL <- function(x,mc=attr(x,"missingClassifier")) {if( is.null(mc) ) is.finite(x)&x<=0.0  else do.call(mc,list(x))=="BDL"  }
is.SZ  <- function(x,mc=attr(x,"missingClassifier")) {if( is.null(mc) ) is.infinite(x)&x<0  else do.call(mc,list(x))=="SZ" }
is.MAR <- function(x,mc=attr(x,"missingClassifier")) {if( is.null(mc) ) is.nan(x)  else do.call(mc,list(x))=="MAR" }
is.MNAR<- function(x,mc=attr(x,"missingClassifier")) {if( is.null(mc) )  is.na(x)&!is.nan(x)  else do.call(mc,list(x))=="MNAR" }
is.NMV <- function(x,mc=attr(x,"missingClassifier")) {if( is.null(mc) )  is.finite(x)&x>0.0 else do.call(mc,list(x))=="NMV" }
is.WZERO <- function(x,mc=attr(x,"missingClassifier")) {is.BDL(x)|is.SZ(x) }
is.WMNAR <- function(x,mc=attr(x,"missingClassifier")) {is.BDL(x)|is.MNAR(x) }

has.missings <- function(x,...) UseMethod("has.missings")

has.missings.default     <- function(x,mc=attr(x,"missingClassifier"),...) {
  (!is.null(x))&&!all(is.NMV(x,mc))
}

has.missings.rmult     <- function(x,mc=attr(x,"missingClassifier"),...) {
  if( is.null(x) )
    return( FALSE )
  if( !is.null(attr(x,"missingProjector")) )
    return( TRUE )
  if( !is.null(orig <- attr(x,"orig")) )
    return( has.missings(orig) )
  else
    !all(is.finite(x))
}
     
getDetectionlimit <- function(x,dl=attr(x,"detectionlimit")) {
  if( is.null(dl) ) dl <- NA
  ifelse(is.finite(x)&x<0,-x,dl)
}

missingType <- function(x,...,mc=attr(x,"missingClassifier"),values=c("NMV","BDT","MAR","MNAR","SZ","Err")) {
  #      finite, infinit, nan, na, x>0&!na  q=(x>0&!na)|nan
  # NA     -       -       -   +    -           -     
  # NaN    -       -       +   +    -           +
  # -Inf   -       +       -   -    -           -         
  # Inf    -       +       -   -    +           +
  # 0/-1   +       -       -   -    -           -
  # +1     +       -       -   -    +           +
  # 6 Faelle -> 3 Bit noetig
  #     infinit, na , !q
  # +1    0      0    0
  # 0/-1  0      0    1
  # NaN   0      1    0
  # NA    0      1    1
  # Inf   1      0    0
  # -Inf  1      0    1
  if( is.null(mc) || identical(mc,missingType)) {
    na <- is.na(x)
    values <- values[c(1,2,3,4,6,5)]
    structure(values[2+is.infinite(x)*4+na*2-((!na&x>0)|is.nan(x))],dim=dim(x),dimnames=dimnames(x))
  }
  else
    do.call(mc,list(x,...))
}



missingSummary <- function(x,...,vlabs=colnames(x),mc=attr(x,"missingClassifier"),values=eval(formals(missingType)$values)) {
  if( is.data.frame(x) )
     x <- data.matrix(x)
  missingType <- structure(c(missingType(x,...,values=1:6)),levels=as.character(values),
            class=c("ordered","factor"))
  # meaning, only faster:
  # missingType <- factor(missingType(x,values=values),levels=values)
  if( length(dim(x)==2 )) {
    if(length(vlabs)!=ncol(x))
      vlabs <- paste("V",1:ncol(x),sep="")
    variable <- structure(rep(1:ncol(x),each=nrow(x)),
                          levels=as.character(vlabs),class="factor")

    as.missingSummary(table(variable,missingType))
  } else {
    as.missingSummary(t(table(missingType)))
  }
}

as.missingSummary <- function(x,...) {
  class(x) <- c("missingSummary",class(x))
  x
}

plot.missingSummary <-function(x,...,main="Missings",legend.text=TRUE,
                    col=c("gray","lightgray","yellow","red","white","magenta")) {
  barplot(t(x),main=main,legend.text=legend.text,...,
          col=col)
}

simulateMissings <- function(x,detectionlimit=NULL,knownlimit=FALSE,MARprob=0.0,MNARprob=0.0,mnarity=0.5,SZprob=0.0) {
  if( is.data.frame(x) )
    x <- data.matrix(x)
  at <- attributes(x)
  x <- unclass(x)
  n <- length(x)
  SZ   <- runif(n)<SZprob
  MAR  <- runif(n)<MARprob
  normal <- function(x) {
    x<- rank(ifelse(is.finite(x),x,1));
    qnorm(x/(length(x)+1))
  }
  w1 <- sqrt(1-mnarity)
  w2 <- sqrt(mnarity)
  MNAR <- pnorm(rnorm(n)*w1+normal(x)*w2) < MNARprob
  x <- clo(x,detectionlimit=detectionlimit,total=NA)
  if( !knownlimit )
    x[is.BDL(x)]<-BDLvalue
  n <- length(x)
  x[MNAR]<-MNARvalue
  x[MAR]<- MARvalue
  x[SZ]<-SZvalue
  attributes(x) <- at
  x
}

zeroreplace <- function(x,d=NULL,a=2/3){
  # ensure objects are either equivalent matrices, data frames or vectors
  W = oneOrDataset(x)
  if(!is.null(d)){
    # if the detection limit is explicitly given, give it the same shape
    Losts = oneOrDataset(is.WZERO(as.matrix(W)),W) # find wide zeroes
    d = oneOrDataset(d,W)
  }else{
    # if not, extract from the data set
    Losts = oneOrDataset(is.BDL(as.matrix(W)),W) # find values below detection limit (encoded as negative)
    d = getDetectionlimit(W)
  }
  # scale down the detection limit...
  d = a*d
  # ... and replace
  W[Losts]=d[Losts]
  attr(W,"Losts")=Losts
  return(W)
}

#zeroreplace <- function(x,d,a=2/3){
#  # ensure all objects are either equivalent matrices, data frames or vectors
#  W = oneOrDataset(x)
#  Losts = oneOrDataset(is.WZERO(as.matrix(W)),W) # find zeroes
#  d = oneOrDataset(d,W)
#  # scale down the detection limit
#  d = a*d
#  # replace
#  W[Losts]=d[Losts]
#  attr(W,"Losts")=Losts
#  return(W)
#}


gsi.mapin01 <- function(i,min=0,max=1) {c(min,min+(max-min)/i,max)}
gsi.mapfrom01 <- function(x) {(x[3]-x[1])/(x[2]-x[1])}
gsi.mapmin <- function(x) {x[1]}
gsi.mapmax <- function(x) {x[3]}

gsi.plots <- list()
gsi.coorInfo <- list()
gsi.setCoorInfo <- function(...) {
  par()
  gsi.coorInfo[[dev.cur()]] <<- list(...)
}
gsi.getCoorInfo <- function() {
  if( dev.cur() <= length(gsi.coorInfo))
    gsi.coorInfo[[dev.cur()]]
  else
    NULL
}

gsi.call  <- function(fkt,...) {
  if( is.character(fkt) )
    do.call(fkt,list(...))
  else
    fkt(...)
}

gsi.add2pairs <- function(x,panel,...,noplot=FALSE) {
  if( dev.cur() <= length(gsi.plots) )
    curplot <- gsi.plots[[dev.cur()]]
  else
    curplot <- NULL
  if( is.null(curplot) ) {
    panel(x[,1],x[,2],...)
  } else {
    if( !missing(panel) )
      curplot$add <- c(curplot$add,list(list(x=x,panel=panel,args=list(...))))
    if(!noplot) do.call("gsi.pairs",curplot)
  }
}




gsi.pairs <- function (x, labels, panel = points, ..., main = NULL, oma = NULL, 
    font.main = par("font.main"), cex.main = par("cex.main"), 
    lower.panel = panel, upper.panel = panel, diag.panel = NULL, 
    text.panel = textPanel, label.pos = 0.5 + has.diag/3, cex.labels = NULL, 
    font.labels = 1, row1attop = TRUE, gap = 1,add=list(),xlim=apply(x,2,range),ylim=apply(x,2,range),log="",noplot=FALSE) 
{
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
        y, txt, cex = cex, font = font)
    localAxis <- function(side, xpd, bg, ...) axis(side, xpd = NA, 
        ...)
    if (!is.matrix(x)) 
        x <- data.matrix(x)
    if (!is.numeric(x)) 
        stop("non-numeric argument to pairs")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
        lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
        upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
        diag.panel <- match.fun(diag.panel)
    if (row1attop) {
        tmp <- lower.panel
        lower.panel <- upper.panel
        upper.panel <- tmp
        tmp <- has.lower
        has.lower <- has.upper
        has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2) 
        stop("only one column in the argument to gsi.pairs")
    has.labs <- TRUE
    if (missing(labels)) {
        labels <- colnames(x)
        if (is.null(labels)) 
            labels <- paste("var", 1:nc)
    }
    else if (is.null(labels)) 
        has.labs <- FALSE
    if( length(dim(xlim)) < 2 ) xlim <- matrix(rep(xlim,ncol(x)),nrow=2)  
    if( length(dim(ylim)) < 2 ) ylim <- matrix(rep(ylim,ncol(x)),nrow=2)  
    if (is.null(oma)) {
        oma <- c(4, 4, 4, 4)
        if (!is.null(main)) 
            oma[3] <- 6
    }

    mycall <- list(x=x,labels=labels,panel=panel,...,main=main,oma=oma,
                   font.main=font.main,cex.main=cex.main,
                   lower.panel=lower.panel, upper.panel=upper.panel,
                   diag.panel=diag.panel,text.panel=text.panel,
                   label.pos=label.pos,cex.labels=cex.labels,
                   font.labels=font.labels,row1attop=row1attop,gap=gap,add=add,
                   xlim=xlim,ylim=ylim,log=log)
    if( noplot ) {
      gsi.plots[[dev.cur()]] <<- mycall
      return(invisible(NULL))
    }
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    for (i in if (row1attop) 
        1:nc
    else nc:1) for (j in 1:nc) {
        plot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
            type = "n", ...,log=ifelse(i==j,"",log),xlim=xlim[,j],ylim=ylim[,i])
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
            box()
            if (i == 1 && (!(j%%2) || !has.upper || !has.lower)) 
                localAxis(1 + 2 * row1attop, ...)
            if (i == nc && (j%%2 || !has.upper || !has.lower)) 
                localAxis(3 - 2 * row1attop, ...)
            if (j == 1 && (!(i%%2) || !has.upper || !has.lower)) 
                localAxis(2, ...)
            if (j == nc && (i%%2 || !has.upper || !has.lower)) 
                localAxis(4, ...)
            mfg <- par("mfg")
            if (i == j) {
                if (has.diag) 
                  diag.panel(as.vector(x[, i]))
                if (has.labs) {
                  par(usr = c(0, 1, 0, 1))
                  if (is.null(cex.labels)) {
                    l.wid <- strwidth(labels, "user")
                    cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                  }
                  text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                    font = font.labels)
                }
            }
            else if (i < j) 
                lower.panel(as.vector(x[, j]), as.vector(x[, 
                  i]), ...)
            else upper.panel(as.vector(x[, j]), as.vector(x[, 
                i]), ...)
            if( i!=j ) for( p in add ) {
              arg <- c(list(p$panel,as.vector(p$x[,j]), as.vector(p$x[,i])),p$args)
              do.call("gsi.call",arg)
            }
            if (any(par("mfg") != mfg)) 
                stop("The panel function made a new plot")
        }
        else par(new = FALSE)
    }
    if (!is.null(main)) 
        mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    gsi.plots[[dev.cur()]] <<- mycall
    invisible(NULL)
}


plot.acomp <- function(x,...,labels=colnames(X),cn=colnames(X),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),margin="acomp",add=FALSE,triangle=!add,col=par("col")) {
  col <- unclass(col)
  X <- oneOrDataset(x)
  oX <- X
  s60 <- sin(pi/3)
  c60 <- cos(pi/3)
  if( NCOL(X) > 3 ) {
    if( margin=="rcomp" )
      infkt <- function(x,y,...) {
        plot.acomp(rcompmargin(X,d=c(gsi.mapfrom01(x),gsi.mapfrom01(y)),pos=1)[,c(3,2,1)],...,aspanel=TRUE,center=center,scale=scale,col=col)
      }
    else if(margin=="acomp") {
      infkt <- function(x,y,...) {
        plot.acomp(acompmargin(X,d=c(gsi.mapfrom01(x),gsi.mapfrom01(y)),pos=1)[,c(3,2,1)],...,aspanel=TRUE,center=center,scale=scale,col=col)
      }
      
    } else {
      if( !is.numeric(margin))
        margin <- match(margin,colnames(X))
      fest <- X[,margin,drop=FALSE]
      X    <- X[,-margin]
      infkt <- function(x,y,...) {
        plot.acomp(acomp(cbind(X[,c(gsi.mapfrom01(y),gsi.mapfrom01(x))],fest)),...,aspanel=TRUE,center=center,scale=scale,col=col)
      }
    }
    nn <- NCOL(X)
    if( add )
      gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
    else
      gsi.pairs(sapply(1:NCOL(X),gsi.mapin01),labels=labels,panel=infkt,...)
  } else {
    if( is.null(cn) ) {
      cn <- c("x","y","z")
    }
    
    
    if( aspanel ) {
      usr <- par("usr"); on.exit(par(usr))
      par( usr=c(0,1,0,1),pty="s" )
      lines(x=c(0,c60,1,0),y=c(0,s60,0,0))
      text(0,0.2,cn[1],pos=4,offset=0.01,xpd=TRUE)
      text(1,0.2,cn[2],pos=2,offset=0.01,xpd=TRUE)
      text(0.5,s60,cn[3],pos=3,offset=0.01,xpd=TRUE)
    } else {
      if( !add ) {
        usr <- par("pty"); on.exit(par(usr))
        par( pty="s" )
        plot(x=c(0,c60,1,0),y=c(0,s60,0,0),
             xlim=c(0,1),ylim=c(0,1),type="n",xlab="",ylab="",
             axes=FALSE)
        gsi.plots[[dev.cur()]]<<-NULL
      }
      if( triangle ) {
        segments(x0=c(0,1,c60),y0=c(0,0,s60),x1=c(1,c60,0),y1=c(0,s60,0))
        #axis(1,pos=0)
        mtext(cn[1],side=1,adj=0,line=1.5)
        mtext(cn[2],side=1,adj=1,line=1.5)
        text(0.5,s60*1.05,cn[3],pos=3,offset=0.01,xpd=TRUE)
      }
    }
    X <- acomp(X,c(1,2,3))
    Y <- scale.acomp(X,center=center,scale=scale)
    gsi.setCoorInfo(mean=if(center) -mean(acomp(X)) else acomp(c(1,1,1)),
                    scale=if(scale) 1/msd(X) else 1)
    x <- Y[,2]+Y[,3]*c60
    y <- Y[,3]*s60
    points(x,y,...,col=col)
    if( id ) {
      if( is.null(idlabs) )
        idlabs <- paste(cn[1],"=",round(X[,1],2),",\n",
                        cn[2],"=",round(X[,2],2),",\n",
                        cn[3],"=",round(X[,3],2))
      return( identify(x,y,idlabs,col=idcol,xpd=NA))
    }
  }
  if( pca && !aspanel ) {
    pca.d <- acomp(princomp(acomp(oX))$Loadings[1,])
    pca.c <- mean(acomp(oX))
    if( center ) pca.c <- pca.c - pca.c
    straight.acomp(pca.c,pca.d,col=col.pca)
  }
  
  return( invisible(NULL))
}

isoPortionLines <- function(...) {
  UseMethod("isoPortionLines",gsi.getCoorInfo()$mean)
}

isoProportionLines <- function(...) {
  UseMethod("isoProportionLines",gsi.getCoorInfo()$mean)
}



isoPortionLines.acomp <- function(by=0.2,at=seq(0,1,by=by),...,parts=1:3,total=1,labs=TRUE,lines=TRUE,unit="") {
  coor <- gsi.getCoorInfo()
  if( coor$scale != 1 )
    stop("Scaling not implemented  in isoPortionLines.acomp")
  for(k in parts) {
    isoCollaps <- function(kkw,k) {
      s <- sum(kkw[-k])
      rcomp(switch(k,c(kkw[k],0,s),c(s,kkw[k],0),c(0,s,kkw[k])))
    }
    dir   <- rep(0,3)
    dir[k]<- 0
    dir[-k]<-c(1,-1)
    dir = rmult(rcomp(dir*unclass(coor$mean))) # directions must be also perturbed!
    for(p in at) {
      start <- rep((1-p)/2,3)
      start[k] <- p
      tx <- rep(0,3)    # to place the text, it is easy to...
      tx[k] <- p         # ... take points on the sides of the diagram
      tx[c(3,1,2)[k]] <- 1-p
      if( p>0 && p<1 ) {
        kw <- clo(start*unclass(coor$mean)) # perturbation must be positive
        tx <- clo(tx*unclass(coor$mean))
        if(lines) straight.rcomp(kw,dir,...)
        kw <- isoCollaps(tx,k)
        if( labs ){
            text(kw[2]+cos(60*pi/180)*kw[3],kw[3]*sin(60*pi/180),
            paste(p*total,unit[(k-1)%%length(unit)+1]),pos=c(2,1,4)[k],...,xpd=TRUE)
           }
      }
    }
  }
  invisible(NULL)
}

#isoPortionLines.acomp <- function(by=0.2,at=seq(0,1,by=by),...,parts=1:3,total=1,labs=TRUE,lines=TRUE,unit="") {
#  coor <- gsi.getCoorInfo()
#  if( coor$scale != 1 )
#    stop("Scaling not implemented  in isoPortionLines.acomp")
#  for(k in parts) {
#    isoCollaps <- function(kw,k) {
#      s <- sum(kw[-k])
#      rcomp(switch(k,c(kw[k],0,s),c(s,kw[k],0),c(0,s,kw[k])))
#    }
#    dir   <- rep(0,3)
#    dir[k]<- 0
#    dir[-k]<-c(1,-1)
#    for(p in at) {
#      start <- rep((1-p)/2,3)
#      start[k] <- p
#      if( p>0 && p<1 ) {
#        kw<-rcomp(acomp(start)-coor$mean)
#        if(lines) straight(kw,rmult(dir),...)
#        kw <- isoCollaps(kw,k)
#        if( labs ) text(kw[2]+cos(60*pi/180)*kw[3],kw[3]*sin(60*pi/180),
#             paste(p*total,unit[(k-1)%%length(unit)+1]),...,pos=c(2,1,4)[k],xpd=TRUE)
#      }
#    }
#  }
#  invisible(NULL)
#}

isoProportionLines.acomp <- function(by=0.2,at=seq(0,1,by=by),...,parts=1:3,labs=TRUE,lines=TRUE) {
  coor <- gsi.getCoorInfo()
  for(k in parts) {
    dir   <- rep(0.25,3)
    dir[k]<- 0.5
    dir = acomp(dir) # * coor$scale # not needed
    for(p in at) {
      if( p>0 && p<1) {
        start <- acomp(switch(k,c(1,1-p,p),c(p,1,1-p),c(1-p,p,1)))
#        start <- acomp(switch(k,c(1E-17,1-p,p),c(p,1E-17,1-p),c(1-p,p,1E-17)))
        kw<-(start+coor$mean)*(coor$scale)  # perturbation must be positive
        if(lines) straight(kw,dir,...)
        kw[k]<- 1E-17
        kw <- acomp(kw)
        if( labs ) text(kw[2]+cos(60*pi/180)*kw[3],kw[3]*sin(60*pi/180),
             paste(p),...,pos=c(4,2,1)[k],xpd=TRUE)
      }
    }
  }
  invisible(NULL)
}

#isoProportionLines.acomp <- function(by=0.2,at=seq(0,1,by=by),...,parts=1:3,labs=TRUE,lines=TRUE) {
#  coor <- gsi.getCoorInfo()
#  for(k in parts) {
#    dir   <- rep(0.25,3)
#    dir[k]<- 0.5
#    for(p in at) {
#      if( p>0 && p<1) {
#        start <- acomp(switch(k,c(1,1-p,p),c(p,1,1-p),c(1-p,p,1)))
#        kw<-(start-coor$mean)*coor$scale
#        if(lines) straight(kw,acomp(dir),...)
#        kw[k]<- 1E-17
#        kw <- acomp(kw)
#        if( labs ) text(kw[2]+cos(60*pi/180)*kw[3],kw[3]*sin(60*pi/180),
#             paste(p),...,pos=c(4,2,1)[k],xpd=TRUE)
#      }
#    }
#  }
#  invisible(NULL)
#}


isoPortionLines.rcomp <- function(by=0.2,at=seq(0,1,by=by),...,parts=1:3,total=1,labs=TRUE,lines=TRUE,unit="") {
  coor <- gsi.getCoorInfo()
  if( coor$scale != 1 || norm(acomp(coor$mean))>0.001 )
    stop("Scaling and centering not implemented  in isoPortionLines.rcomp")
  for(k in parts) {
    isoCollaps <- function(kw,k) {
      s <- sum(kw[-k])
      rcomp(switch(k,c(kw[k],0,s),c(s,kw[k],0),c(0,s,kw[k])))
    }
    dir   <- rep(0,3)
    dir[k]<- 0
    dir[-k]<-c(1,-1)
    for(p in at) {
      start <- rep((1-p)/2,3)
      start[k] <- p
      if( p>0 && p<1 ) {
        try({
        kw <- rcomp(start) # isoCollaps(rcomp(start)-coor$mean)
        if(lines ) straight(kw,rmult(dir),...)
        kw <- isoCollaps(kw,k)
        if( labs ) text(kw[2]+cos(60*pi/180)*kw[3],kw[3]*sin(60*pi/180),
             paste(p*total,unit[(k-1)%%length(unit)+1]),...,pos=c(2,1,4)[k],xpd=TRUE)
      },silent=FALSE)
      }
    }
  }
  invisible(NULL)
}


isoProportionLines.rcomp <- function(by=0.2,at=seq(0,1,by=by),...,parts=1:3,labs=TRUE,lines=TRUE) {
  coor <- gsi.getCoorInfo()
  if( coor$scale != 1 || norm(acomp(coor$mean))>0.001)
    stop("Scaling not implemented  in isoPortionLines.rcomp")
  for(k in parts) {
    dir   <- rep(0.25,3)
    dir[k]<- 0.5
    for(p in at) {
      if( p>0 && p<1) {
        start <- acomp(switch(k,c(1,1-p,p),c(p,1,1-p),c(1-p,p,1)))
        kw<-start
        if(lines) straight(kw,acomp(dir),...)
        kw[k]<- 1E-17
        kw <- acomp(kw)
        if( labs ) text(kw[2]+cos(60*pi/180)*kw[3],kw[3]*sin(60*pi/180),
             paste(p),...,pos=c(4,2,1)[k],xpd=TRUE)
      }
    }
  }
  invisible(NULL)
}



plot.rcomp <- function(x,...,labels=colnames(X),cn=colnames(X),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),margin="rcomp",add=FALSE,col=par("col")) {
if( center || scale )
  warning("Scaling and centring meaningless for rcomp-compositions");
X <- oneOrDataset(x)
oX<-X
s60 <- sin(pi/3)
c60 <- cos(pi/3)
if( NCOL(X) > 3 ) {
  if( margin=="rcomp" )
    infkt <- function(x,y,...) {
      plot.rcomp(rcompmargin(X,d=c(gsi.mapfrom01(x),gsi.mapfrom01(y)),pos=1)[,c(3,2,1)],...,aspanel=TRUE,col=col)
    }
  else if(margin=="acomp") {
    infkt <- function(x,y,...) {
      plot.rcomp(acompmargin(X,d=c(gsi.mapfrom01(x),gsi.mapfrom01(y)),pos=1)[,c(3,2,1)],...,aspanel=TRUE,col=col)
    }
    
  } else {
    if( !is.numeric(margin))
      margin <- match(margin,colnames(X))
    fest <- X[,margin,drop=FALSE]
    X    <- X[,-margin]
    infkt <- function(x,y,...) {
      plot.rcomp(acomp(cbind(X[,c(gsi.mapfrom01(y),gsi.mapfrom01(x))],fest)),...,aspanel=TRUE)
    }
  }
  if( add )
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  else
    gsi.pairs(sapply(1:NCOL(X),gsi.mapin01),labels=labels,panel=infkt,...)
  } else {
    if( is.null(cn) ) {
      cn <- c("x","y","z")
    }
    X <- rcomp(X,c(1,2,3))
    Y <- X
    gsi.setCoorInfo(mean=rcomp(c(1,1,1)),
                    scale=1)
    x <- Y[,2]+Y[,3]*c60
    y <- Y[,3]*s60
    
    if( aspanel ) {
      usr <- par("usr"); on.exit(par(usr))
      par( usr=c(0,1,0,1),pty="s" )
      lines(x=c(0,c60,1,0),y=c(0,s60,0,0))
      text(0,0.2,cn[1],pos=4,offset=0.01,xpd=TRUE)
      text(1,0.2,cn[2],pos=2,offset=0.01,xpd=TRUE)
      text(0.5,s60,cn[3],pos=3,offset=0.01,xpd=TRUE)
      
    } else {
      if( !add ) {
        usr <- par("pty"); on.exit(par(usr))
        par( pty="s" )
        plot(x=c(0,c60,1,0),y=c(0,s60,0,0),
             xlim=c(min(c(0,x)),max(c(1,x))),
             ylim=c(min(c(0,y)),max(c(1,y))),
             type="n",xlab="",ylab="",
             axes=FALSE)
        gsi.plots[[dev.cur()]]<<-NULL
        
        arrows(x0=c(0,1,c60),y0=c(0,0,s60),x1=c(1,c60,0),y1=c(0,s60,0),angle=15)
        # axis(1,pos=0)
        mtext(cn[1],side=1,adj=0)
        mtext(cn[2],side=1,adj=1)
        text(0.5,s60,cn[3],pos=3,offset=0.01,xpd=TRUE)
      }
    }
    points(x,y,...,col=col)
    if( id ) {
      if( is.null(idlabs) )
        idlabs <- paste(cn[1],"=",round(X[,1],2),",\n",
                        cn[2],"=",round(X[,2],2),",\n",
                        cn[3],"=",round(X[,3],2))
      return( identify(x,y,idlabs,col=idcol,xpd=NA))
    }
  }
if( pca && ! aspanel) {
  pca.d <- princomp(cpt(oX))$loadings[,1]
  pca.c <- mean(rcomp(oX))
  straight.rcomp(pca.c,pca.d,col=col.pca)
}
return( invisible(NULL))
}

plot.aplus <- function(x,...,labels=colnames(X),cn=colnames(X),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),add=FALSE,logscale=TRUE,col=par("col")) {
  col <- unclass(col)
X <- oneOrDataset(x)
oX <- X
if( NCOL(X) > 2 ) {
    infkt <- function(x,y,...) {
      plot.aplus(X[,c(x[1],y[1])],...,aspanel=TRUE,center=center,scale=scale,logscale=logscale,add=add,col=col)
    }
    usr <- par(c("xlog","ylog")); on.exit(par(usr))
    #if( !add ) par(xlog=logscale,ylog=logscale)
    if( add ) 
      gsi.add2pairs(matrix(1:NCOL(X),nrow=1),panel=infkt,...)
    else 
      gsi.pairs(matrix(1:NCOL(X),nrow=1),labels=labels,panel=infkt,...,log=ifelse(logscale,"xy",""))
  } else {
    if( is.null(cn) ) {
      cn <- c("x","y")
    }
    x <- X[,1]
    y <- X[,2]
    if( ! add ) {
      if( aspanel ) {
        usr <- par("usr"); on.exit(par(usr))
        if( logscale )
          par( xlog=TRUE,ylog=TRUE,usr=c(min(log10(x[is.NMV(x)])),max(log10(x[is.NMV(x)])),
                                         min(log10(y[is.NMV(y)])),max(log10(y[is.NMV(y)]))))
        else
          par( usr=c(min(x,na.rm=TRUE),max(x,na.rm=TRUE),min(y,na.rm=TRUE),max(y,na.rm=TRUE)) )
                                        #axis(1)
                                        #axis(2)
      } else {
        if( !add ) {
          plot(x=c(1),y=c(1),
               xlim=range(x),ylim=range(y),type="n",
               log=ifelse(logscale,"xy",""),xlab=cn[1],ylab=cn[2])
          gsi.plots[[dev.cur()]]<<-NULL
        }
      }
    }
    points(x,y,...,col=col)
    if( id ) {
      if( is.null(idlabs) )
        idlabs <- paste(cn[1],"=",round(X[,1],2),",\n",
                        cn[2],"=",round(X[,2],2))
      return( identify(x,y,idlabs,col=idcol,xpd=NA))
    }
  }
if( pca && ! aspanel) {
  pca.d <- ilt.inv(princomp(ilt(oX))$loadings[,1])
  pca.c <- mean(aplus(oX))
  straight.aplus(pca.c,pca.d,col=col.pca)
}
return( invisible(NULL))
}

plot.rplus <- function(x,...,labels=colnames(X),cn=colnames(X),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),add=FALSE,logscale=FALSE,xlim=apply(X,2,function(x) c(0,max(x,na.rm=TRUE))),ylim=xlim,col=par("col")) {
  col<- unclass(col)
if( scale )
  warning("Centering meaningless for rplus-amounts");
X <- oneOrDataset(x)
oX <-X
if( NCOL(X) > 2 ) {
    infkt <- function(x,y,...) {
      plot.rplus(X[,c(x[1],y[1])],...,aspanel=TRUE,center=center,scale=scale,logscale=logscale,col=col)
    }
    if( add )
      gsi.add2pairs(matrix(1:NCOL(X),nrow=1),infkt,...)
    else {
      gsi.pairs(matrix(1:NCOL(X),nrow=1),labels=labels,panel=infkt,...,xlim=xlim,ylim=ylim)
    }
    
  } else {
    if( is.null(cn) ) {
      cn <- c("x","y")
    }
    x <- X[,1]
    y <- X[,2]
    if( aspanel && ! add ) {
      usr <- par("usr"); on.exit(par(usr))
 #     if( logscale )
 #       par( xlog=TRUE,ylog=TRUE,usr=c(log10(min(x)),log10(max(x)),log10(min(y)),log10(max(y))))
 #     else
 #       par( usr=c(0,max(x),0,max(y)) )
 #                                       #axis(1)
 #                                       #axis(2)
    } else {
      if( !add ) {
        plot(x=c(1),y=c(1),
             xlim=xlim[,1],ylim=ylim[,2],type="n",
             log=ifelse(logscale,"xy",""),xlab=cn[1],ylab=cn[2])
        gsi.plots[[dev.cur()]]<<-NULL
      }
    }
    points(x,y,...,col=col)
    if( id ) {
      if( is.null(idlabs) )
        idlabs <- paste(cn[1],"=",round(X[,1],2),",\n",
                        cn[2],"=",round(X[,2],2))
      return( identify(x,y,idlabs,col=idcol,xpd=NA))
    }
  }
if( pca && ! aspanel ) {
  pca.d <- princomp(iit(oX))$loadings[,1]
  pca.c <- mean(rplus(oX))
  straight.rplus(pca.c,pca.d,col=col.pca)
}

return( invisible(NULL))
}

plot.rmult <- function(x,...,labels=colnames(X),cn=colnames(X),aspanel=FALSE,id=FALSE,idlabs=NULL,idcol=2,center=FALSE,scale=FALSE,pca=FALSE,col.pca=par("col"),add=FALSE,logscale=FALSE,col=par("col")) {
X <- oneOrDataset(x)
oX <- X
if( NCOL(X) > 2 ) {
    infkt <- function(x,y,...) {
      plot.rmult(X[,c(x[1],y[1])],...,aspanel=TRUE,center=center,scale=scale,pca=pca,col.pca=col.pca,logscale=logscale,col=col)
    }
    if( add )
      gsi.add2pairs(matrix(1:NCOL(X),nrow=1),infkt,...)
    else
      gsi.pairs(matrix(1:NCOL(X),nrow=1),labels=labels,panel=infkt,...)
  } else {
    if( is.null(cn) ) {
      cn <- c("x","y")
    }
    x <- X[,1]
    y <- X[,2]
    if( aspanel && ! add ) {
      usr <- par("usr"); on.exit(par(usr))
      if( logscale )
        par( xlog=TRUE,ylog=TRUE,usr=c(log10(min(x)),log10(max(x)),log10(min(y)),log10(max(y))))
      else
        par( usr=c(min(x),max(x),min(y),max(y)) )
                                        #axis(1)
                                        #axis(2)
    } else {
      if( !add ) {
        plot(x=c(1),y=c(1),
             xlim=range(x),ylim=range(y),type="n",
             log=ifelse(logscale,"xy",""),xlab=cn[1],ylab=cn[2])
        gsi.plots[[dev.cur()]]<<-NULL
      }
    }
    points(x,y,...,col=col)
    if( id ) {
      if( is.null(idlabs) )
        idlabs <- paste(cn[1],"=",round(X[,1],2),",\n",
                        cn[2],"=",round(X[,2],2))
      return( identify(x,y,idlabs,col=idcol,xpd=NA))
    }
  }
if( pca && ! aspanel ) {
  pca.d <- iit.inv(princomp(iit(oX))$loadings[,1])
  pca.c <- mean(oX)
  straight.rmult(pca.c,pca.d,col=col.pca)
}
return( invisible(NULL))
}


lines.acomp <- function(x,...,steps=30) {
  X <- oneOrDataset(x)
  if( ncol(X) > 3 ) {
    infkt <- function(x,y,...) {
      
      lines.acomp(gsi.plotmargin(X,c(gsi.mapfrom01(x),gsi.mapfrom01(y)),
                                 get("margin",environment(gsi.plots[[dev.cur()]]$panel))),
                                 ...,steps=steps)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    Y <- X[-1,,drop=FALSE]
    X <- X[-nrow(X),,drop=FALSE]
    s60 <- sin(pi/3)
    c60 <- cos(pi/3)
    l   <- rep((0:steps)/steps,NROW(X))
    i   <- rep(1:NROW(X),each=steps+1)
    XP  <- unclass(ilr.inv((1-l)*ilr(X[i,,drop=FALSE]) +
                           l*ilr(Y[i,,drop=FALSE])))
    x <- XP[,2]+XP[,3]*c60
    y <- XP[,3]*s60
    lines(x,y,...)
  }
}

lines.rcomp <- function(x,...,steps=30) {
  X <- oneOrDataset(x)
  if( ncol(X) > 3 ) {
    infkt <- function(x,y,...) {
      
      lines.rcomp(gsi.plotmargin(X,c(gsi.mapfrom01(x),gsi.mapfrom01(y)),
                                 get("margin",environment(gsi.plots[[dev.cur()]]$panel))),
                  ...,steps=steps)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    Y <- X[-1,,drop=FALSE]
    X <- X[-nrow(X),,drop=FALSE]
    s60 <- sin(pi/3)
    c60 <- cos(pi/3)
    l   <- rep(c((0:steps)/steps),NROW(X))
    i   <- rep(1:NROW(X),each=steps+1)
    XP  <- unclass(convex.rcomp(X[i,,drop=FALSE],Y[i,,drop=FALSE],l))
    x <- XP[,2]+XP[,3]*c60
    y <- XP[,3]*s60
    lines(x,y,...)
  }
}


lines.aplus <- function(x,...,steps=30) {
  X <- oneOrDataset(x)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      lines.aplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],...)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    Y <- X[-1,,drop=FALSE]
    X <- X[-nrow(X),,drop=FALSE]
    l   <- rep((0:steps)/steps,NROW(X))
    i   <- rep(1:NROW(X),each=steps+1)
    XP  <- unclass(ilt.inv((1-l)*ilt(X[i,,drop=FALSE]) + l*ilt(Y[i,,drop=FALSE])))
    x <- XP[,1]
    y <- XP[,2]
    lines(x,y,...)
  }
}

lines.rplus <- function(x,...,steps=30) {
  X <- oneOrDataset(x)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      lines.rplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],...)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    Y <- X[-1,,drop=FALSE]
    X <- X[-nrow(X),,drop=FALSE]
    l   <- rep((0:steps)/steps,NROW(X))
    i   <- rep(1:NROW(X),each=steps+1)
    XP  <- unclass(iit.inv((1-l)*iit(X[i,,drop=FALSE]) + l*iit(Y[i,,drop=FALSE])))
    x <- XP[,1]
    y <- XP[,2]
    lines(x,y,...)
  }
}

lines.rmult <- function(x,...,steps=30) {
  X <- oneOrDataset(x)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      lines.rmult(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],...)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    Y <- X[-1,,drop=FALSE]
    X <- X[-nrow(X),,drop=FALSE]
    l   <- rep((0:steps)/steps,NROW(X))
    i   <- rep(1:NROW(X),each=steps+1)
    XP  <- (1-l)*unclass(X)[i,,drop=FALSE] + l*unclass(Y)[i,,drop=FALSE]
    x <- XP[,1]
    y <- XP[,2]
    lines(x,y,...)
  }
}


segments.acomp <- function(x0,y,...,steps=30) {
  X <- oneOrDataset(x0,y)
  Y <- oneOrDataset(y,x0)
  if( ncol(X) > 3 ) {
    infkt <- function(x,y,...) {
      mar <- get("margin",environment(gsi.plots[[dev.cur()]]$panel))
      segments.acomp(
                     gsi.plotmargin(X,c(gsi.mapfrom01(x),gsi.mapfrom01(y)),mar),
                     gsi.plotmargin(Y,c(gsi.mapfrom01(x),gsi.mapfrom01(y)),mar)
                                 ,
                                 ...,steps=steps)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    s60 <- sin(pi/3)
    c60 <- cos(pi/3)
    l   <- rep(c((0:steps)/steps,NA),NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    XP  <- unclass(ilr.inv((1-l)*ilr(X[i,]) + l*ilr(Y[i,])))
    x <- XP[,2]+XP[,3]*c60
    y <- XP[,3]*s60
    lines(x,y,...)
  }
}

segments.rcomp <- function(x0,y,...,steps=30) {
  X <- oneOrDataset(x0,y)
  Y <- oneOrDataset(y,x0)
  if( ncol(X) > 3 ) {
    infkt <- function(x,y,...) {
      mar <- get("margin",environment(gsi.plots[[dev.cur()]]$panel))
      segments.rcomp(
                     gsi.plotmargin(X,c(gsi.mapfrom01(x),gsi.mapfrom01(y)),mar),
                     gsi.plotmargin(Y,c(gsi.mapfrom01(x),gsi.mapfrom01(y)),mar)
                                 ,
                                 ...,steps=steps)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  }
  s60 <- sin(pi/3)
  c60 <- cos(pi/3)
  l   <- rep(c((0:steps)/steps,NA),NROW(X))
  i   <- rep(1:NROW(X),each=steps+2)
  XP  <- unclass(convex.rcomp(X[i,],Y[i,],l))
  x <- XP[,2]+XP[,3]*c60
  y <- XP[,3]*s60
  lines(x,y,...)
}

segments.aplus <- function(x0,y,...,steps=30) {
  X <- oneOrDataset(x0,y)
  Y <- oneOrDataset(y,x0)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      segments.aplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],
                     Y[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],...)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    l   <- rep(c((0:steps)/steps,NA),NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    XP  <- unclass(ilt.inv((1-l)*ilt(X[i,]) + l*ilt(Y[i,])))
    x <- XP[,1]
    y <- XP[,2]
    lines(x,y,...)
  }
}

segments.rplus <- function(x0,y,...,steps=30) {
  X <- oneOrDataset(x0,y)
  Y <- oneOrDataset(y,x0)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      segments.rplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],
                     Y[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],...)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    l   <- rep(c((0:steps)/steps,NA),NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    XP  <- unclass(iit.inv((1-l)*iit(X[i,]) + l*iit(Y[i,])))
    x <- XP[,1]
    y <- XP[,2]
    lines(x,y,...)
  }
}

segments.rmult <- function(x0,y,...,steps=30) {
  X <- oneOrDataset(x0,y)
  Y <- oneOrDataset(y,x0)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      segments.rmult(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],
                     Y[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))],...)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    l   <- rep(c((0:steps)/steps,NA),NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    XP  <- (1-l)*unclass(X)[i,] + l*unclass(Y)[i,]
    x <- XP[,1]
    y <- XP[,2]
    lines(x,y,...)
  }
}


#segments.panel.rcomp <- function(X,Y,...,steps=30) {
#  function(what,...) {
#    if( !is.null(colnames(what)) )
#      what <- colnames(what)
#    X <- gsi.margin(X,what)
#    Y <- gsi.margin(Y,what)
#    segments.rcomp(X,Y,...,steps=steps)
#  }
#}

#segments.panel.aplus <- function(X,Y,...,steps=30) {
#  function(what,...) {
#    if( !is.null(colnames(what)) )
#      what <- colnames(what)
#    X <- X[,what]
#    Y <- Y[,what]
#    segments.aplus(X,Y,...,steps=steps)
#  }
#}

#segments.panel.rplus <- function(X,Y,...,steps=30) {
#  function(what,...) {
#    if( !is.null(colnames(what)) )
#      what <- colnames(what)
#    X <- X[,what]
#    Y <- Y[,what]
#    segments.rplus(X,Y,...,steps=steps)
#  }
#}

#segments.panel.rmult <- function(X,Y,...,steps=30) {
#  function(what,...) {
#    if( !is.null(colnames(what)) )
#      what <- colnames(what)
#    X <- X[,what]
#    Y <- Y[,what]
#    segments.rplus(X,Y,...,steps=steps)
#  }
#}



gsi.closespread <- function(spread) {
  if(length(dim(spread))>3) {
    return(apply(spread,1,gsi.closespread))
  }
  d <- nrow(spread)
  Pmat <- diag(d)-1/d
  row.names(Pmat)<-row.names(spread)
  Pmat %*% spread %*% t(Pmat)
}

gsi.spreadToIsoSpace <- function(spread) {
  if(length(dim(spread))>3) {
    return(apply(spread,1,gsi.closespread))
  }
  d <- nrow(spread)
  V <- ilrBase(D=d)
  t(V) %*% spread %*% V
}

ellipses <- function(mean,...) UseMethod("ellipses",mean)

ellipses.rcomp <- function(mean,var,r=1,...,steps=360) {
mean <- oneOrDataset(ipt(mean))
sp <- var
w  <- seq(0,2*pi,length.out=steps)
for(i in 1:nrow(mean)) {
    if( length(dim(var))==3 )
      sp<-var[i,,]
    isp <- gsi.spreadToIsoSpace(sp)
    mi <- mean[i,]
    eisp <- eigen(isp,TRUE)
    X <- t(mi+ t(sqrt(eisp$values[1])*r*cos(w) %o% eisp$vectors[,1] +
             sqrt(eisp$values[2])*r*sin(w) %o% eisp$vectors[,2])) 
    lines.rcomp(ucipt.inv(X),...)
  }
}

ellipses.rplus <- function(mean,var,r=1,...,steps=360) {
mean <- oneOrDataset(iit(mean))
sp <- var
w  <- seq(0,2*pi,length.out=steps)
for(i in 1:nrow(mean)) {
    if( length(dim(var))==3 )
      sp<-var[i,,]
#    isp <- gsi.spreadToIsoSpace(sp)
    mi <- mean[i,]
    eisp <- eigen(sp,TRUE)
    X <- t(mi+ t(sqrt(eisp$values[1])*r*cos(w) %o% eisp$vectors[,1] +
             sqrt(eisp$values[2])*r*sin(w) %o% eisp$vectors[,2])) 
    lines.rmult(X,...)
  }
}

ellipses.rmult <- function(mean,var,r=1,...,steps=360) {
mean <- oneOrDataset(mean)
sp <- var
w  <- seq(0,2*pi,length.out=steps)
for(i in 1:nrow(mean)) {
  if( length(dim(var))==3 )
    sp<-var[i,,]
  mi <- mean[i,]
  eisp <- eigen(sp,TRUE)
  X <- t(mi+ t(sqrt(eisp$values[1])*r*cos(w) %o% eisp$vectors[,1] +
               sqrt(eisp$values[2])*r*sin(w) %o% eisp$vectors[,2])) 
  lines.rmult(X,...)
  }
}

ellipses.acomp <- function(mean,var,r=1,...,steps=360) {
mean <- ilr(mean)
sp <- var
w  <- seq(0,2*pi,length.out=steps)
for(i in 1:nrow(mean)) {
    if( length(dim(var))==3 )
      sp<-var[i,,]
    isp <- gsi.spreadToIsoSpace(sp)
    mi <- unclass(mean)[i,]
    eisp <- eigen(isp,TRUE)
    X <- t(mi + t(sqrt(eisp$values[1])*r*cos(w) %o% eisp$vectors[,1] +
             sqrt(eisp$values[2])*r*sin(w) %o% eisp$vectors[,2])) 
    lines.acomp(ilr.inv(X),...)
  }
}

ellipses.aplus <- function(mean,var,r=1,...,steps=360) {
mean <- oneOrDataset(ilt(mean))
sp <- var
w  <- seq(0,2*pi,length.out=steps)
for(i in 1:nrow(mean)) {
    if( length(dim(var))==3 )
      sp<-var[i,,]
   # isp <- gsi.spreadToIsoSpace(sp)
    mi <- unclass(mean)[i,]
    eisp <- eigen(sp,TRUE)
    X <- t(mi + t(sqrt(eisp$values[1])*r*cos(w) %o% eisp$vectors[,1] +
             sqrt(eisp$values[2])*r*sin(w) %o% eisp$vectors[,2])) 
    lines.aplus(ilt.inv(X),...)
  }
}

straight  <- function(x,...) UseMethod("straight",x)

straight.acomp <- function(x,d,...,steps=30) {
  X <- oneOrDataset(x,d)
  d <- oneOrDataset(d,x)
  if( ncol(X) > 3 ) {
    infkt <- function(x,y,...) {
      mar <- get("margin",environment(gsi.plots[[dev.cur()]]$panel))
      straight.acomp(
                     gsi.plotmargin(acomp(X),c(gsi.mapfrom01(x),gsi.mapfrom01(y)),mar),
                     gsi.plotmargin(acomp(d),c(gsi.mapfrom01(x),gsi.mapfrom01(y)),mar)
                                 ,
                                 ...,steps=steps)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    d <- normalize(acomp(d)) 
    X <- perturbe(X,power.acomp(d,-scalar(acomp(X),acomp(d)))) 
    s60 <- sin(pi/3)
    c60 <- cos(pi/3)
    l   <- rep(2*c((0:steps)/steps,NA)-1,NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    XP  <- acomp(clr.inv(clr(X[i,]) + (2*l)^3*clr(d[i,])))
    x <- XP[,2]+XP[,3]*c60
    y <- XP[,3]*s60
    lines(x,y,...)
  }
}

straight.aplus <- function(x,d,...,steps=30) {
  X <- oneOrDataset(x,d)
  d <- oneOrDataset(d,x)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      straight.aplus(
                     aplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))]),
                     aplus(d[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))]),
                                 ...,steps=steps)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {
    s <- log(d[,2])/log(d[,1])
#  if( par("xlog") && par("ylog") )
#    abline(exp(log(x[,2])-log(x[,1])/s),s,untf=FALSE,...)
#  else {
    l   <- rep(2*c((0:steps)/steps,NA)-1,NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    r   <- par("usr")
    if( ! par("xlog") ) r[1:2] <- log(c(r[2]/100,r[2]))
    if( ! par("ylog") ) r[3:4] <- log(c(r[4]/100,r[4]))
    r   <- abs(r[2]-r[1])+abs(r[4]-r[3])
    XP  <- aplus(X[i,]) + r*l*normalize(aplus(d[i,]))
    lines.rmult(XP,...) # lines.aplus(XP,...) produced double lines
                                        #    warning("straight.aplus not yet implemented in nonlog coordinates");
                                        #  }
  }
}


straight.rcomp <- function(x,d,...,steps=30) {
  X <- oneOrDataset(x,d)
  d <- oneOrDataset(d,x)
  l1 <- apply(-X/d,1,function(x) {max(x[x<=0])})*0.99999 
  l2 <- apply(-X/d,1,function(x) {min(x[x>=0])})*0.99999
  X1 <- rcomp(rmult(X)+(l1*rmult(d))) 
  X2 <- rcomp(rmult(X)+(l2*rmult(d)))
  if( ncol(X) > 3 ) {
    infkt <- function(x,y,...) {
      mar <- get("margin",environment(gsi.plots[[dev.cur()]]$panel))
      segments.rcomp(
                     gsi.plotmargin(rcomp(X1),c(gsi.mapfrom01(x),gsi.mapfrom01(y)),mar),
                     gsi.plotmargin(rcomp(X2),c(gsi.mapfrom01(x),gsi.mapfrom01(y)),mar)
                     ,
                     ...,steps=steps)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else { 
    segments.rcomp(X1,X2,...)
  }
}
#straight.rplus <- function(x,d,...,steps=30) {
#  x <- oneOrDataset(x,d)
#  d <- oneOrDataset(d,x)
#  s <- d[,2]/d[,1]
#  abline(x[,2]-x[,1]/s,s,untf=TRUE)
#  l1 <- apply(-x/d,1,function(x) {max(c(-100,x[x<=0]))}) 
#  l2 <- apply(-x/d,1,function(x) {min(c(100,x[x>=0]))})
#  X1 <- rcomp(gsi.add(x,gsi.mul(l1,d))) 
#  X2 <- rcomp(gsi.add(x,gsi.mul(l2,d)))
#  segments.rplus(X1,X2,...)
#}

straight.rplus <- function(x,d,...,steps=30) {
  X <- oneOrDataset(x,d)
  d <- oneOrDataset(d,x)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      straight.rplus(
                     rplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))]),
                     rmult(d[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))]),
                                 ...,steps=steps)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {  
    #s <- log(d[,2])/log(d[,1])
                                        #  if( par("xlog") && par("ylog") )
                                        #    abline(exp(log(x[,2])-log(x[,1])/s),s,untf=FALSE,...)
                                        #  else {
    l   <- rep(2*c((0:steps)/steps,NA)-1,NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    r   <- par("usr")
    if( par("xlog") ) r[1:2] <- 10^(r[1:2])
    if( par("ylog") ) r[3:4] <- 10^(r[3:4])
    r   <- abs(r[2]-r[1])+abs(r[4]-r[3])
    XP  <- rmult(X[i,]) + r*l*normalize(rmult(d[i,]))
    lines.rmult(XP,...)
  }
}

straight.rmult <- function(x,d,...,steps=30) {
  X <- oneOrDataset(x,d)
  d <- oneOrDataset(d,x)
  if( ncol(X) > 2 ) {
    infkt <- function(x,y,...) {
      straight.rmult(
                     aplus(X[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))]),
                     aplus(d[,c(gsi.mapfrom01(x),gsi.mapfrom01(y))]),
                                 ...,steps=steps)
    }
    gsi.add2pairs(sapply(1:NCOL(X),gsi.mapin01),infkt,...)
  } else {  
    #s <- log(d[,2])/log(d[,1])
                                        #  if( par("xlog") && par("ylog") )
                                        #    abline(exp(log(x[,2])-log(x[,1])/s),s,untf=FALSE,...)
#  else {
    l   <- rep(2*c((0:steps)/steps,NA)-1,NROW(X))
    i   <- rep(1:NROW(X),each=steps+2)
    r   <- par("usr")
    if( par("xlog") ) r[1:2] <- 10^(r[1:2])
    if( par("ylog") ) r[3:4] <- 10^(r[3:4])
    r   <- abs(r[2]-r[1])+abs(r[4]-r[3])
    XP  <- rmult(X[i,]) + r*l*normalize(rmult(d[i,]))
    lines.rmult(XP,...)
                                        #    warning("straight.aplus not yet implemented in nonlog coordinates");
                                        #  }
  }
}


#straight.panel.acomp <- function(x,d,...,steps=30) {
#  function(what,...) {
#    if( !is.null(colnames(what)) )
#      what <- colnames(what)
#    x <- margin(x,what)
#    d <- margin(d,what)
#    straight.acomp(x,d,...,steps=steps)
#  }
#}

#straight.panel.rcomp <- function(x,d,...,steps=30) {
#  function(what,...) {
#    if( !is.null(colnames(what)) )
#      what <- colnames(what)
#    x <- margin(x,what)
#    d <- margin(d,what)
#    straight.rcomp(x,d,...,steps=steps)
#  }
#}



rDirichlet.acomp <- function(n,alpha) {
  acomp(sapply(alpha,rgamma,n=n))
}

rDirichlet.rcomp <- function(n,alpha) {
  rcomp(sapply(alpha,rgamma,n=n))
}



runif.acomp <- function(n,D) rDirichlet.acomp(n,rep(1,D))
runif.rcomp <- function(n,D) rDirichlet.rcomp(n,rep(1,D))

rnorm.aplus <- function(n,mean,var) {
  D <- NCOL(oneOrDataset(mean))
  perturbe.aplus(ilt.inv(matrix(rnorm(n*length(mean)),ncol=D) %*% chol(var)),
                mean)
}

dnorm.aplus <- function(x,mean,var) {
  x <- aplus(x)
  mean <- aplus(mean)
  w <- ilt(x-mean)
  D <- ncol(oneOrDataset(x))
  if( length(dim(w)) == 2 ) 
    u <- c(rep(1,ncol(w))%*%t((solve(var,t(w)))*t(w)))
  else
    u <- sum(solve(var,w)*w)
  exp(-u/2)/sqrt(2^D*pi^D*det(var))
}

dnorm.rmult <- function(x,mean,var) {
  w <- rmult(x)-rmult(mean)
  D <- gsi.getD(x)
  if( length(dim(w)) == 2 ) 
    u <- c(rep(1,ncol(w))%*%t((solve(var,t(w)))*t(w)))
  else
    u <- sum(solve(var,w)*w)
  exp(-u/2)/sqrt(2^D*pi^D*det(var))
}

rnorm.acomp <- function(n,mean,var) {
  D <- NCOL(oneOrDataset(mean))
  perturbe(ilr.inv(matrix(rnorm(n*(D-1)),ncol=D-1) %*%
                         chol(clrvar2ilr(var))),
                mean)
}

dnorm.acomp <- function(x,mean,var) {
  x <- acomp(x)
  mean <- acomp(mean)
  w <- ilr(x-mean)
  D <- ncol(oneOrDataset(x))
  if( length(dim(w)) == 2 ) 
    u <- c(rep(1,D-1)%*%t((solve(clrvar2ilr(var),t(w)))*t(w)))
  else
    u <- sum(solve(clrvar2ilr(var),w)*w)
  exp(-u/2)/sqrt(2*pi*det(clrvar2ilr(var)))
}



rlnorm.rplus <- function(n,meanlog,varlog) {
  D <- NCOL(oneOrDataset(meanlog))
  rplus(perturbe.aplus(exp(matrix(rnorm(n*length(meanlog)),ncol=D) %*% chol(varlog)),
                exp(meanlog)))
}

dlnorm.rplus <- function(x,meanlog,varlog) {
  xx <- oneOrDataset(x)
  w <- ilt(x)-meanlog
  if( length(dim(w)) == 2 ) {
    u <- c(rep(1,ncol(w))%*%t((solve(varlog,t(w)))*t(w)))
    v <- c(exp(log(xx) %*% rep(1,ncol(xx)))) 
  }
  else {
    u <- solve(varlog,w)%*%w
    v <- prod(x)
  }
  exp(-u/2)/sqrt(2*pi*det(varlog))/v
}


rnorm.rplus <- function(n,mean,var) {
  D <- ncol(var)
  rplus(pmax(rmult(matrix(rnorm(n*D),ncol=D) %*% chol(var))+rmult(mean),0))
}

rnorm.rmult <- function(n,mean,var) {
  D <- ncol(var)
  rmult(matrix(rnorm(n*D),ncol=D) %*% chol(var))+rmult(mean)
}

rnorm.rcomp <- function(n,mean,var) {
  D <- ncol(var)
  rcomp(pmax(ilr2clr(matrix(rnorm(n*D),ncol=D-1) %*% chol(clrvar2ilr(var)))+rplus(rcomp(mean)),0))
}



gsi.margin <- function(X,...) UseMethod("gsi.margin")

gsi.plotmargin <- function(X,d,margin) {
  X <- oneOrDataset(X)
  if( margin=="rcomp" )
    rcompmargin(X,d=d,pos=1)[,c(3,2,1)]
  else if( margin=="acomp" )
    acompmargin(X,d=d,pos=1)[,c(3,2,1)]
  else {
    if( ! is.numeric(margin))
      margin <- match(margin,colnames(X))
    fest <- X[,margin,drop=FALSE]
    X    <- X[,-margin]
    acomp(cbind(fest,X[,d]))[,c(3,2,1)]
  }
}
  
gsi.margin.acomp <- function(X,what,...,margin="acomp") {
  if( margin == "sub" )
    acomp(X,what)
  else if( margin=="rcomp" )
      rcompmargin(X,what)
  else if( margin=="acomp")
    acompmargin(X,what)
  else {
    if( !is.numeric(what) )
      what <- match(what,colnames(X))
    if( !is.numeric(margin))
      margin <- match(margin,colnames(X))
    acomp(X,c(what,margin))
  }
}

gsi.margin.rcomp <- function(X,what,...,margin="rcomp") {
  if( margin == "sub" )
    acomp(X,what)
  else if( margin=="rcomp" )
    rcompmargin(X,what)
  else if( margin=="acomp")
    acompmargin(X,what)
  else {
    if( !is.numeric(what) )
      what <- match(what,colnames(X))
    if( !is.numeric(margin))
      margin <- match(margin,colnames(X))
    rcomp(X,c(what,margin))
  }
}

gsi.margin.aplus <- function(X,what,...) {
  aplus(X,what)
}

gsi.margin.rplus <- function(X,what,...) {
  rplus(X,what)
}

gsi.isSingleRow <- function(X) {
  return( NROW(X) == 1 || NCOL(X) ==1 )
}

barplot.acomp <- function(height,...,legend.text=TRUE,beside=FALSE,total=1) {
  X <- height
  if( gsi.isSingleRow(X) )
    barplot(t(rbind(gsi.plain(acomp(X,total=total)),0)),c(1,0),...,legend.text=legend.text,
            beside=beside)
  else
    barplot(gsi.plain(t(acomp(X,total=total))),...,legend.text=legend.text,beside=beside)
}

barplot.rcomp <- function(height,...,legend.text=TRUE,beside=FALSE,total=1) {
  X <- height
  if( gsi.isSingleRow(X) )
    barplot(t(rbind(gsi.plain(rcomp(X,total=total)),0)),c(1,0),...,legend.text=legend.text,beside=beside)
  else
    barplot(gsi.plain(t(rcomp(X,total=total))),...,legend.text=legend.text,beside=beside);
}

barplot.aplus <- function(height,...,legend.text=TRUE,beside=TRUE) {
  X <- height
  if( gsi.isSingleRow(X) )
    barplot(t(rbind(gsi.plain(aplus(X)),0)),c(1,0),...,legend.text=legend.text,
            beside=beside)
  else
    barplot(gsi.plain(t(aplus(X))),...,legend.text=legend.text,beside=beside);
}

barplot.rplus <- function(height,...,legend.text=TRUE,beside=TRUE) {
  X <- height
  if( gsi.isSingleRow(X) )
    barplot(t(rbind(gsi.plain(rplus(X)),0)),c(1,0),...,
            legend.text=legend.text,beside=beside)
  else
    barplot(gsi.plain(t(rplus(X))),...,legend.text=legend.text,beside=beside);
}



split.acomp <- function(x,f,drop=FALSE,...) {
  cls <- class(x)
  lapply(split(1:NROW(x),f,drop=drop,...),function(i) structure(x[i,],class=cls))
}
split.rcomp <- split.acomp
split.aplus <- split.acomp
split.rplus <- split.acomp
split.rmult <- split.acomp

as.data.frame.acomp <- function(x,...) as.data.frame.matrix(unclass(x))
as.data.frame.rcomp <- function(x,...) as.data.frame.matrix(unclass(x))
as.data.frame.aplus <- function(x,...) as.data.frame.matrix(unclass(x))
as.data.frame.rplus <- function(x,...) as.data.frame.matrix(unclass(x))
as.data.frame.rmult <- function(x,...) as.data.frame.matrix(unclass(x))

gsi.addclass <- function(x,cls) {
  class(x) <- c(cls,attr(x,"class"))
  x
}

gsi <- function() {}


princomp.acomp <- function(x,...,scores=TRUE) {
  cl <- match.call()
  D <- gsi.getD(x)
  tmp <- princomp(clr(x),...,scores=scores)
  tmp$sdev        <- tmp$sdev[-D]
  tmp$loadings    <- structure(tmp$loadings[,-D],class="loadings")
  tmp$Center      <- clr.inv(tmp$center)
  tmp$Loadings  <- acomp(clr.inv(t(tmp$loadings)),total=D)
  tmp$DownLoadings<- acomp(clr.inv(t(-tmp$loadings)),total=D)
  tmp$call <- cl
  gsi.addclass(tmp,"princomp.acomp")
}

print.princomp.acomp <- function(x,...) {
  NextMethod("print",x,...)
  cat("Mean (compositional):\n")
  print(x$Center)
  cat("+Loadings (compositional):\n")
  print(x$Loadings)
  cat("-Loadings (compositional):\n")
  print(x$DownLoadings)
  invisible(x)
}

plot.princomp.acomp <- function(x,y=NULL,...,npcs=min(10,length(x$sdev)),
                                        type=c("screeplot","variance",
                                          "biplot","loadings","relative"),
                                main=NULL,scale.sdev=1) {
  if( missing(main) ) main <- deparse(substitute(x))
  type <- match.arg(type)
  if( type=="biplot" )
    biplot(x,...,main=main)
  else if( type=="loadings" ) {
    if( is.na(scale.sdev) )
      scl<-1
    else
      scl<-scale.sdev*x$sdev
    barplot(acomp((x$Loadings*scl)[1:npcs,]),...,
            main=main,total=gsi.getD(x$Loadings))
  }
  else if(type=="relative") {
    tmp <- relativeLoadings(x,scale.sdev=scale.sdev)[,1:npcs]
    tmp <- barplot(t(tmp),...,main=main,beside=TRUE,legend=TRUE)
    abline(h=1)
    invisible(tmp)
  } else  {
  if( type=="screeplot" ) type <- "lines"
  if( type=="variance" ) type <- "barplot"
  screeplot(x,...,npcs=npcs, main=main,type=type)
  }
}

predict.princomp.acomp <- function(object,newdata,...) {
  NextMethod("predict",object,newdata=clr(newdata),...)
}

                                
#panel.princomp.acomp <- function(x,choice,t,...){
#  straight.panel.acomp(x$Center,x$Loadings)
#}


princomp.rcomp <- function(x,...,scores=TRUE) {
  cl <- match.call()
  D <- gsi.getD(x)
  tmp <- princomp(cpt(x),...,scores=scores)
  tmp$sdev        <- tmp$sdev[-D]
  tmp$loadings    <- structure(tmp$loadings[,-D],class="loadings")
  tmp$Center      <- cpt.inv(tmp$center)
  tmp$Loadings    <- rmult(t(tmp$loadings))
  tmp$call <- cl
  gsi.addclass(tmp,"princomp.rcomp")
}


print.princomp.rcomp <- function(x,...) {
  NextMethod("print",x,...)
}

plot.princomp.rcomp <- function(x,y=NULL,...,npcs=min(10,length(x$sdev)),
                                        type=c("screeplot","variance",
                                          "biplot","loadings","relative"),
                                main=NULL,scale.sdev=1) {
  if( missing(main) ) main <- deparse(substitute(x))
  type <- match.arg(type)
  if( type=="biplot" )
    biplot(x,...,main=main)
  else if( type=="loadings" ) {
    if( is.na(scale.sdev) )
      scl<-1
    else
      scl<-scale.sdev*x$sdev
    barplot((rmult(t(x$loadings))*scl)[1:npcs,],...,
            main=main)
  }
  else if(type=="relative") {
    tmp <- relativeLoadings(x,scale.sdev=scale.sdev)[,1:npcs]
    tmp <- barplot(t(tmp),...,main=main,beside=TRUE,legend=TRUE)
    invisible(tmp)
  } else  {
  if( type=="screeplot" ) type <- "lines"
  if( type=="variance" ) type <- "barplot"
  screeplot(x,...,npcs=npcs, main=main,type=type)
  }
}

predict.princomp.rcomp <- function(object,newdata,...) {
  NextMethod("predict",object,newdata=cpt(newdata),...)
}

princomp.aplus <- function(x,...,scores=TRUE) {
  cl <- match.call()
  D <- gsi.getD(x)
  tmp <- princomp(ilt(x),...,scores=scores)
  tmp$Center      <- ilt.inv(tmp$center)
  tmp$Loadings  <- ilt.inv(t(tmp$loadings))
  tmp$DownLoadings<- ilt.inv(t(-tmp$loadings))
  tmp$call <- cl
  gsi.addclass(tmp,"princomp.aplus")
}

print.princomp.aplus <- function(x,...) {
  NextMethod("print",x,...)
  cat("Mean (compositional):\n")
  print(x$Center)
  cat("+Loadings (compositional):\n")
  print(x$Loadings)
  cat("-Loadings (compositional):\n")
  print(x$DownLoadings)
  invisible(x)
}

plot.princomp.aplus <- function(x,y=NULL,...,npcs=min(10,length(x$sdev)),
                                        type=c("screeplot","variance",
                                          "biplot","loadings","relative"),
                                main=NULL,scale.sdev=1) {
  if( missing(main) ) main <- deparse(substitute(x))
  type <- match.arg(type)
  if( type=="biplot" )
    biplot(x,...,main=main)
  else if( type=="loadings" ) {
        if( is.na(scale.sdev) )
      scl<-1
    else
      scl<-scale.sdev*x$sdev
    barplot(aplus((x$Loadings*scl)[1:npcs,]),...,
            main=main)
  }
  else if(type=="relative") {
    tmp <- relativeLoadings(x,scale.sdev=scale.sdev)[,1:npcs]
    tmp <- barplot(t(tmp),...,main=main,beside=TRUE,legend=TRUE)
    abline(h=1)
    invisible(tmp)
  } else  {
  if( type=="screeplot" ) type <- "lines"
  if( type=="variance" ) type <- "barplot"
  screeplot(x,...,npcs=npcs, main=main,type=type)
  }
}

predict.princomp.aplus <- function(object,newdata,...) {
  NextMethod("predict",object,newdata=ilt(newdata),...)
}

princomp.rplus <- function(x,...,scores=TRUE) {
  cl <- match.call()
  tmp <- princomp(iit(x),...,scores=scores)
  tmp$Center      <- iit.inv(tmp$center)
  tmp$call <- cl
  tmp$Loadings <- rmult(t(tmp$loadings))
  gsi.addclass(tmp,"princomp.rplus")
}

print.princomp.rplus <- function(x,...) {
  NextMethod("print",x,...)
  cat("Mean:\n")
  print(x$Center)
  cat("Loadings:\n")
  print(x$Loadings)
  invisible(x)
}

plot.princomp.rplus <- function(x,y=NULL,...,npcs=min(10,length(x$sdev)),
                                        type=c("screeplot","variance",
                                          "biplot","loadings","relative"),
                                main=NULL,scale.sdev=1) {
  if( missing(main) ) main <- deparse(substitute(x))
  type <- match.arg(type)
  if( type=="biplot" )
    biplot(x,...,main=main)
  else if( type=="loadings" ) {
    if( is.na(scale.sdev) )
      scl<-1
    else
      scl<-scale.sdev*x$sdev
    barplot((rmult(t(x$loadings))*scl)[1:npcs,],...,
            main=main)
  }
  else if(type=="relative") {
    tmp <- relativeLoadings(x,scale.sdev=scale.sdev)[,1:npcs]
    tmp <- barplot(t(tmp),...,main=main,beside=TRUE,legend=TRUE)
    invisible(tmp)
  } else  {
  if( type=="screeplot" ) type <- "lines"
  if( type=="variance" ) type <- "barplot"
  screeplot(x,...,npcs=npcs, main=main,type=type)
  }
}

predict.princomp.rplus <- function(object,newdata,...) {
  NextMethod("predict",object,newdata=iit(newdata),...)
}


princomp.rmult <- function(x,...) {
  princomp(unclass(x),...)
}

gsi.pairrelativeMatrix <- function(names) {
  n  <- length(names)
  ii <- rep(1:n,n)
  jj <- rep(1:n,each=n)
  jgi <- jj>ii
  ii <- ii[jgi]
  jj <- jj[jgi]
  N <- length(ii)
  erg <- rep(0,N*n)
  erg[1:N+N*(ii-1)]<-  1
  erg[1:N+N*(jj-1)]<- -1
  erg <- matrix(erg,nrow=N,ncol=n)
  colnames(erg)  <- names
  row.names(erg) <- paste(names[ii],names[jj],sep="/") 
  erg
}

relativeLoadings <- function(x,...) UseMethod("relativeLoadings",x)
relativeLoadings.princomp.acomp <- function(x,...,log=FALSE,scale.sdev=NA,cutoff=0.1) {
  if( is.na(scale.sdev) ) 
    scale <- rep(1,length(x$sdev))
  else 
    scale <- scale.sdev*x$sdev
  bl <- gsi.pairrelativeMatrix(row.names(x$loadings)) %*% (unclass(x$loadings) %*% diag(scale))
  colnames(bl) <- colnames(x$loadings)
  if( !log )
    bl <- exp(bl)
  structure(bl,class="relativeLoadings.princomp.acomp",cutoff=cutoff,scale=scale,log=log)
}

relativeLoadings.princomp.aplus <- function(x,...,log=FALSE,scale.sdev=NA,cutoff=0.1) {
  if( is.na(scale.sdev) ) 
    scale <- rep(1,length(x$sdev))
  else 
    scale <- scale.sdev*x$sdev
  bl <- gsi.pairrelativeMatrix(row.names(x$loadings)) %*% (unclass(x$loadings) %*% diag(scale))
  colnames(bl) <- colnames(x$loadings)
  if( !log )
    bl <- exp(bl)
  structure(bl,class="relativeLoadings.princomp.aplus",cutoff=cutoff,scale=scale,log=log)
}

relativeLoadings.princomp.rcomp <- function(x,...,scale.sdev=NA,cutoff=0.1) {
  if( is.na(scale.sdev) ) 
    scale <- rep(1,length(x$sdev))
  else 
    scale <- scale.sdev*x$sdev
  bl <- gsi.pairrelativeMatrix(row.names(x$loadings)) %*% (unclass(x$loadings) %*% diag(scale))
  colnames(bl) <- colnames(x$loadings)
  structure(bl,class="relativeLoadings.princomp.rcomp",cutoff=cutoff,scale=scale)
}

relativeLoadings.princomp.rplus <- function(x,...,scale.sdev=NA,cutoff=0.1) {
  if( is.na(scale.sdev) ) 
    scale <- rep(1,length(x$sdev))
  else 
    scale <- scale.sdev*x$sdev
  bl <- gsi.pairrelativeMatrix(row.names(x$loadings)) %*% (unclass(x$loadings) %*% diag(scale))
  colnames(bl) <- colnames(x$loadings)
  structure(bl,class="relativeLoadings.princomp.rplus",cutoff=cutoff,scale=scale)
}


print.relativeLoadings.princomp.acomp <- function(x,...,
                                                 cutoff=attr(x,"cutoff"),
                                                 digits=2
                                                 ) {
  bt <- format(x,digits=digits)
  if( !attr(x,"log") ) { 
    if( !is.na(cutoff) )
      bt[abs(log(x))<cutoff] <- ""
  } else {
    if( !is.na(cutoff) )
      bt[abs(x)<cutoff] <- ""
  }
  print(bt,quote=FALSE)
  x
}

print.relativeLoadings.princomp.aplus <- print.relativeLoadings.princomp.acomp

print.relativeLoadings.princomp.rcomp <- function(x,...,
                                                 cutoff=attr(x,"cutoff"),
                                                 digits=2
                                                 ) {
  bt <- format(x,digits=digits)
  bt[abs(x)<cutoff] <- ""
  print(bt,quote=FALSE)
  x
}

print.relativeLoadings.princomp.rplus <- function(x,...,
                                                 cutoff=attr(x,"cutoff"),
                                                 digits=2
                                                 ) {
  bt <- format(x,digits=digits)
  bt[abs(x)<cutoff] <- ""
  print(bt,quote=FALSE)
  x
}


plot.relativeLoadings.princomp.acomp<- function(x,...) {
  barplot(t(x),...,beside=TRUE)
  if( !attr(x,"log") )
    abline(h=1)
  invisible(x)
}
plot.relativeLoadings.princomp.aplus<- function(x,...) {
  barplot(t(x),...,beside=TRUE)
  if( !attr(x,"log") )
    abline(h=1)
  invisible(x)
}

plot.relativeLoadings.princomp.rcomp<- function(x,...) {
  barplot(t(x),...,beside=TRUE)
  invisible(x)
}

plot.relativeLoadings.princomp.rplus<- function(x,...) {
  barplot(t(x),...,beside=TRUE)
  invisible(x)
}



read.geoeas <- function (file){
#read title
print("reading title:...")
title <- read.table(file, nrows=1, quote="", colClasses="character", sep="\t")
print("reading title: OK")

#read number of variables
print("reading number of variables:...")
nvars <- read.table(file, skip=1, nrows=1, sep="\t")
nvars <- as.integer(nvars)
#read labels of the variables
print("reading variables:...")
labels <- scan(file, skip=2, nmax=nvars, nlines= nvars, sep="\t", quote="", what="character")

#read data matrix
print("reading dataset:...")
dades <- read.table(file, skip=2+nvars)
print("reading dataset: OK")
#relate variables with their names
colnames(dades)=as.vector(labels)
#add title as an attribute
attr(dades,"title") <- as.character(title)
return(dades)
}

read.geoEAS <- function(file){ read.geoeas(file) }



#read.standard <- function (file){
##read title
#print("reading title:...")
#title <- read.table(file, nrows=1, quote="", colClasses="character", sep="\t")
#print("reading title: OK")

#read data matrix
#print("reading dataset:...")
#dades <- read.table(file, skip=1, header=T)
#print("reading dataset: OK")#

#add title as an attribute
#attr(dades,"title") <- as.character(title)
#return(dades)
#}

 #  =========================================================================#
# Based on the the tedrahedron plot from MixeR
# http://vlado.fmf.uni-lj.si/pub/MixeR/
# Vladimir Batagelj & Matevz Bren
# plots mix object 'm' into tetrahedron and exports it in kinemage format
#   clu - partition -> color of points
#   vec - vector of values -> size of points
#   col - color of points if clu=NULL
#   scale - relative size of points
#   king - FALSE (for Mage); TRUE (for King)
kingTetrahedron <- function(X,parts=1:4,file="tmptetrahedron.kin",clu=NULL,vec=NULL,king=TRUE,scale=0.2,col=1,title="Compositional Tetrahedron"){
  m <- list()
  m$mat <- data.frame(clo(X,parts=parts))
  m$tit <- title
  kinColors <- c('white', 'red', 'blue', 'yellow',
    'green', 'cyan', 'magenta', 'pink', 'lime',
    'orange', 'peach', 'gold', 'purple', 'sea',
    'gray', 'sky', 'brown', 'lilac', 'hotpink',
    'yellowtint', 'pinktint', 'peachtint',
    'greentint', 'bluetint', 'lilactint',
    'deadwhite', 'deadblack', 'invisible')
  warn <- ""
  if (king) warn <- "\n*** works only with KiNG viewer: http://kinemage.biochem.duke.edu/"
  head <- paste("@text\n",
    file," / ",date(),
    "\nDataset: ", m$tit,warn,
"\nLayout obtained using composition
based on  MixeR
http://vlado.fmf.uni-lj.si/pub/MixeR/
Vladimir Batagelj & Matevz Bren
Institute of Mathematics, Physics and Mechanics
Ljubljana, Slovenia
@kinemage 1
@caption\n", m$tit,
"\n@fontsizeword 18
@zoom 1.00
@thinline
@zclipoff
@group {Complete} dominant animate movieview = 1
@spherelist 0 radius= 0.20\n",sep="")
  foot <-
"@vectorlist {}  color=  blue
P 9.500, 0.500, 9.500
0.500, 9.500, 9.500
P 9.500, 0.500, 9.500
0.500, 0.500, 0.500
P 9.500, 0.500, 9.500
9.500, 9.500, 0.500
P 0.500, 9.500, 9.500
0.500, 0.500, 0.500
P 0.500, 9.500, 9.500
9.500, 9.500, 0.500
P 0.500, 0.500, 0.500
9.500, 9.500, 0.500\n"
  kin <- file(file,"w")
  n <- nrow(m$mat)
  if (is.null(clu)) {clu <- rep(col,n)}
  clu <- c(clu,0,0,0,0)
  if (is.null(vec)) {vec <- rep(1,n)}
  vec <- c(vec,1,1,1,1)
  size <- scale/max(vec)
  rn <- c(dimnames(m$mat)[[1]],"A","B","C","D")
  a <- 10/9+c(m$mat[,1],1,0,0,0)
  b <- c(m$mat[,2],0,1,0,0)
  c <- c(m$mat[,3],0,0,1,0)
  d <- c(m$mat[,4],0,0,0,1)
  x <- (a - b - c + d)*0.45
  y <- (a - b + c - d)*0.45
  z <- (a + b - c - d)*0.45
  cat(head,file=kin)
  for (i in 1:length(rn)) {
    color <- "white"
    if (clu[i]>0) color <- kinColors[2+(clu[i]-1)%%18]
    cat("{",rn[i],"} ", color," ",file=kin,sep="")
    if (king) cat(" r=", vec[i]*size," ",file=kin,sep=" ")
    cat(10*x[i],10*(1-y[i]),10*z[i],"\n",file=kin,sep=" ")
  }
  cat(foot,file=kin)
  close(kin)
}

