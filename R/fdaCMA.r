### filename: fdaCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 3.10.2007; major update and revisions: 23.10.2007
#
### Brief description:
#
#  Returns an object of class cloutput.
#
### Further comments and notes:
#  - The computation algorithm is based on Ripley (1996)
#  - for computation, CCA is used
#  - see also flexDA, where Generalized Additive Models are used
#    via mgcv.
#
###**************************************************************************###

setGeneric("fdaCMA", function(X, y, f, learnind, comp=1, plot = FALSE)
           standardGeneric("fdaCMA"))


### signature X=matrix, y=numeric, f=missing:

setMethod("fdaCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, comp=1, plot=FALSE){
nrx <- nrow(X)
ly <- length(y)
if(nrx != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(ncol(X) >= length(y)) stop("Too many variables selected. \n")
if(missing(learnind)) learnind <- 1:nrx
if(length(learnind) > nrx)
stop("length of 'learnind' must be smaller than the number of observations. \n")
y <- as.factor(y)
levels(y) <- 1:nlevels(y)
maxcomp <- min(nlevels(y)-1, ncol(X))
if(comp > maxcomp){
 warning("'comp' must be <= min(number of predictors, number of classes-1); set to maximum possible value \n")
 comp <- maxcomp
}
if(nlevels(y) > 2) mode <- "multiclass"
else mode <- "binary"
y <- as.numeric(y)-1
Ylearn <- y[learnind]
G <- model.matrix(~as.factor(Ylearn)-1)
taby <- colSums(G)
Xlearn <- X[learnind,,drop=FALSE]
#Xcenter <- scale(Xlearn, scale=FALSE)
#svdX <- svd(Xcenter)
#S <- sqrt(nrow(Xlearn))*svdX$v %*% diag(1/svdX$d)
#Xscale <- Xcenter %*% S
CCA <- cancor(Xlearn, G)
A <- CCA$xcoef[,1:comp,drop=FALSE]
Xsubs <- Xlearn %*% A
centroids <- scale(t(Xsubs) %*% G, FALSE, taby)
Xtest <- X[-learnind,,drop=FALSE]
if(plot) yy <- y
if(nrow(Xtest) == 0){
y <- Ylearn
distm <- t(apply(Xsubs, 1, function(z) colSums((z - centroids)^2)))
yhat <- apply(distm, 1, which.min)-1
distm <- distm - rowMeans(distm)
prob <- exp(-distm)
prob <- prob/rowSums(prob)
}
else{ y <- y[-learnind]
      #Xtestcenter <- scale(Xtest, scale=FALSE)
      #svdXtest <- svd(Xtestcenter)
      #Stest <- sqrt(length(y))*svdXtest$v %*% diag(1/svdXtest$d)
      #Xscaletest <- Xtestcenter %*% Stest
      Xsubstest <- Xtest %*% A
      #Xsubstest <- Xtestcenter %*% A
      distm <- t(apply(Xsubstest, 1, function(z) colSums((z - centroids)^2)))
      yhat <- apply(distm, 1, which.min)-1
      distm <- distm - rowMeans(distm)
      prob <- exp(-distm)
      prob <- prob/rowSums(prob)
    }
if(plot){
 if(comp > 2) stop("If 'plot=TRUE', number of components have to be 1 or 2. \n")
 else{
    if(exists("Xsubstest")){
    XX <- matrix(nrow=length(yy), ncol=comp)
    XX[learnind,] <- Xsubs
    XX[-learnind,] <- Xsubstest
    }
    else XX <- Xsubs
    pchs <- ifelse(1:length(yy) %in% learnind, 15, 17)
    if(comp == 1){
      plot(XX[,1], rep(0.5, length(yy)), pch=pchs, col=yy+2, xlab="first canonical variate",
           main="learning set: square     test set: triangle    centroids: circles",
           cex.main=1, ylab="")
      points(centroids, rep(0.5, ncol(G)),col=(1:ncol(G))+1, pch=16, cex=2.5)
    }
    else{
      plot(XX[,1], XX[,2], pch=pchs, col=yy+2, xlab="first canonical variate",
           main="learning set: square     test set: triangle    centroids: circles",
           cex.main=1, ylab="second canonical variate")
      points(centroids[1,], centroids[2,], col=(1:ncol(G))+1, pch=16, cex=2.5)
    }
 }
}
new("cloutput", yhat=yhat, y=y, learnind = learnind, prob = prob, method = "FDA", mode=mode)
})

### signature X=matrix, y=factor, f=missing:

setMethod("fdaCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, comp=1, plot=FALSE){
fdaCMA(X, y=as.numeric(y)-1, learnind=learnind, comp=comp, plot=plot)
})

### signature X=data.frame, y=missing f=formula:

setMethod("fdaCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f,learnind, comp=1, plot=FALSE){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
fdaCMA(as.matrix(X), y=as.numeric(y)-1, learnind=learnind,
      comp = comp, plot=plot)  })

### signature X=ExpressionSet, y=character

setMethod("fdaCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind, comp=1, plot=FALSE){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          fdaCMA(X=X, y=y, learnind=learnind, comp = comp, plot = plot)})