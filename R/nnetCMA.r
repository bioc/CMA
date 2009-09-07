### filename: nnetCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski, adapted from A-L Boulesteix
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 18.9.2007
#
### Brief description:
#
#  Returns an object of class cloutput.
#
### Further comments and notes:
#
#
#
###**************************************************************************###

### generic

setGeneric("nnetCMA", function(X, y, f, learnind, eigengenes = FALSE, ...) standardGeneric("nnetCMA"))

### X=matrix, y=numeric, f=missing method

setMethod("nnetCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, eigengenes = FALSE, ...){
require(nnet, quietly=TRUE)
nrx <- nrow(X)
ly <- length(y)
if(nrx != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) learnind <- 1:nrx
if(length(learnind) > nrx)
stop("length of 'learnind' must be smaller than the number of observations. \n")
y <- as.factor(y)
levels(y) <- 1:nlevels(y)
dotsCall <- substitute(list(...))
ll <- eval(dotsCall)
if(!hasArg(size)) ll$size <- 1
if(!hasArg(MaxNWts)) ll$MaxNWts <- 1000
if(!hasArg(trace)) ll$trace <- FALSE
if(nlevels(y) > 2) mode <- "multiclass"
else mode <- "binary"
y <- as.numeric(y)-1
Ylearn <- y[learnind]
G <- class.ind(as.factor(Ylearn))

noweights <- (ncol(X)+1)*(ll$size)+(ll$size+1)*ncol(G)
if(noweights > ll$MaxNWts) stop("Number of weights too large. Either increase
                                 'MaxNwts' (s. package nnet, function nnet) or
                                 perform a variable selection \n")
Xlearn <- X[learnind,,drop = FALSE]

#svd using learning set only
if(eigengenes){
      svdX <- svd(Xlearn)
      svalue <- svdX$d
      svaluePos <- seq(svalue)[svalue > 0]
      svalue <- svalue[svaluePos]
      Xlearn <- svdX$u[, svaluePos] %*% diag(svalue)
}



ll$x <- Xlearn ; ll$y <- G
output.nnet <- do.call("nnet", args = ll)
Xtest <- X[-learnind,,drop=FALSE]

#using training-svd for testdata 
if(eigengenes){
   Xtest<- Xtest %*% svdX$v
   colnames(Xtest)<-1:ncol(Xtest)
}


if(nrow(Xtest) == 0){ Xtest <- Xlearn ; y <- Ylearn }
else y <- y[-learnind]
pred.test <- predict(object=output.nnet, newdata=Xtest)
yhat <- apply(pred.test, 1, which.max)-1
new("cloutput", yhat=yhat, y=y, learnind = learnind,
     prob = pred.test, method = "nnet", mode=mode)
})

### signature X=matrix, y=factor, f=missing:

setMethod("nnetCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, eigengenes = FALSE, ...){
nnetCMA(X, y=as.numeric(y)-1, learnind=learnind, eigengenes = eigengenes, ...)
})

### signature X=data.frame, f=formula

setMethod("nnetCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, eigengenes = FALSE, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
nnetCMA(as.matrix(X), y=y, learnind=learnind, eigengenes = eigengenes, ...)})


### signature: X=ExpressionSet, y=character.

setMethod("nnetCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind, eigengenes = FALSE, ...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          nnetCMA(X=X, y=y, learnind=learnind, eigengenes = eigengenes, ...)})
