### filename: knnCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski, adapted from A-L Boulesteix
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 16.10.2007
#
### Brief description:
#
#   'Usual' nearest neighbour classifier.
#   Returns an object of class cloutput.
#
### Further comments and notes:
#
###**************************************************************************###

setGeneric("knnCMA", function(X, y, f, learnind, models=FALSE, ...)
           standardGeneric("knnCMA"))

### X=matrix, y=numeric, f=missing

setMethod("knnCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, models=FALSE,...){
require(class, quietly=TRUE)
nrx <- nrow(X)
ly <- length(y)
if(nrx != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) learnind <- 1:nrx
if(length(learnind) > nrx)
stop("length of 'learnind' must be smaller than the number of observations. \n")
y <- as.factor(y)
levels(y) <- 1:nlevels(y)
if(nlevels(y) > 2) mode <- "multiclass"
else mode <- "binary"
y <- as.numeric(y)-1
Ylearn <- y[learnind]
Xlearn <- X[learnind,]
Xtest <- X[-learnind,,drop=FALSE]
if(nrow(Xtest) == 0) stop("For k nearest neighbours, a test set must be given \n")
else{
y <- y[-learnind]
yhat <- knn(train = Xlearn, test = Xtest, cl=Ylearn, ...)
}


if(models==TRUE)
	modd<-list(NULL)
if(models==FALSE)
	modd<-list(NULL)

new("cloutput", yhat=as.numeric(yhat)-1, y=y, learnind = learnind,
     prob = matrix(data = NA, nrow = length(learnind), ncol=length(unique(yhat))),
     method = "knn", mode=mode,model=modd)})

### signature X=matrix, y=factor

setMethod("knnCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, models=FALSE,...){
knnCMA(X, y=as.numeric(y)-1, learnind=learnind,models=models,...)
})

### signature X=data.frame, f=formula

setMethod("knnCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, models=FALSE,...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
knnCMA(as.matrix(X), y=y, learnind=learnind,models=models,...)})


### signature: X=ExpressionSet, y=character.

setMethod("knnCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind,models=FALSE,...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          knnCMA(X=X, y=y, learnind=learnind, models=models,...)})