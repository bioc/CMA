### filename: gbmCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski, adapted from A-L Boulesteix
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 16.10.2007
#
### Brief description:
#
#  Returns an object of class cloutput.
#
### Further comments and notes:
#
#   - Boosting based on (low-order) decision trees
#   - Different losses can be used.
#
#
###**************************************************************************###

### generic

setGeneric("gbmCMA", function(X, y, f, learnind, ...)
           standardGeneric("gbmCMA"))

### X=matrix, y=numeric, f=missing.

setMethod("gbmCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, ...){
require(gbm, quietly=TRUE)
nrx <- nrow(X)
ly <- length(y)
if(nrx != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) learnind <- 1:nrx
if(length(learnind) > nrx)
stop("length of 'learnind' must be smaller than the number of observations. \n")
y <- as.factor(y)
levels(y) <- 1:nlevels(y)
if(nlevels(y) > 2) stop("'gbmCMA' only possible for binary classification \n")
else mode <- "binary"
y <- as.numeric(y)-1
dotsCall <- substitute(list(...))
ll <- eval(dotsCall)
if(!hasArg(n.minobsinnode))ll$n.minobsinnode <- 1
if(!hasArg(bag.fraction)) ll$bag.fraction <- 1
if(!hasArg(n.trees)) ll$n.trees <- 1000
if(!hasArg(verbose)) ll$verbose <- FALSE
ll$x <- X[learnind,,drop=FALSE]
ll$y <- y[learnind]
output.gbm <- do.call(gbm.fit, args=ll)
Xtest <- data.frame(X[-learnind,,drop=FALSE])
if(nrow(Xtest) == 0){ Xtest <- X[learnind,,drop=FALSE] ; y <- y[learnind] }
else y <- y[-learnind]
prob <- predict(output.gbm, newdata=Xtest, n.trees=ll$n.trees, type="response")
prob <- cbind(1-prob, prob)
yhat <- apply(prob, 1, which.max)-1
new("cloutput", yhat=yhat, y=y, learnind = learnind,
     prob = prob, method = "gbm", mode=mode)
})

#### signature X=matrix, y=numeric, f=missing

setMethod("gbmCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, ...){
gbmCMA(X, y=as.numeric(y)-1, learnind=learnind,...)
})

### signature X=data.frame, f=formula

setMethod("gbmCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
gbmCMA(as.matrix(X), y=y, learnind=learnind,...)})


### signature: X=ExpressionSet, y=character.

setMethod("gbmCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind,...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          gbmCMA(X=X, y=y, learnind=learnind, ...)})