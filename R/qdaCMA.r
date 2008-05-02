### filename: qdaCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski
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


setGeneric("qdaCMA", function(X, y, f, learnind, ...)
           standardGeneric("qdaCMA"))


### signature X=matrix, y=numeric, f=missing:

setMethod("qdaCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, ...){
require(MASS, quietly=TRUE)
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
if(nlevels(y) > 2) mode <- "multiclass"
else mode <- "binary"
y <- as.numeric(y)-1
Ylearn <- y[learnind]
Xlearn <- data.frame(X[learnind,])
output.qda <- qda(Xlearn, grouping=factor(Ylearn),...)
Xtest <- data.frame(X[-learnind,,drop=FALSE])
if(nrow(Xtest) == 0){ Xtest <- Xlearn ; y <- Ylearn }
else y <- y[-learnind]
pred.test <- predict(object=output.qda, newdata=data.frame(Xtest))
new("cloutput", yhat=as.numeric(pred.test$class)-1, y=y, learnind = learnind,
     prob = pred.test$posterior, method = "QDA", mode=mode)
})

### signature X=matrix, y=factor, f=missing:

setMethod("qdaCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, ...){
qdaCMA(X, y=as.numeric(y)-1, learnind=learnind,...)
})

### signature X=data.frame, f=formula

setMethod("qdaCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
qdaCMA(as.matrix(X), y=y, learnind=learnind,...)})


### signature: X=ExpressionSet, y=character.

setMethod("qdaCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind,...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          qdaCMA(X=X, y=y, learnind=learnind, ...)})