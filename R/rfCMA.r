### Classification based on random forests (using the package randomForest)
###
### This function uses the functions randomForest and predict.randomForest from
### the package randomForest by Liaw and Wiener
###
### filename: rfCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski, adapted from A-L Boulesteix
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 26.9.2007
#
### Brief description:
#
#  Returns an object of class cloutput.
#
### Further comments and notes:
#
###**************************************************************************###

### generic

setGeneric("rfCMA", function(X, y, f, learnind, varimp = TRUE, seed = 111, models=FALSE,...)
           standardGeneric("rfCMA"))

### X=matrix, y=numeric, f=missing - method

setMethod("rfCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, varimp = TRUE, seed=111, models=FALSE,...){
require(randomForest, quietly=TRUE)
nrx <- nrow(X)
ly <- length(y)
if(nrx != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) learnind <- 1:nrx
y <- as.factor(y)
if(nlevels(y) > 2) mode <- "multiclass"
else mode <- "binary"
levels(y) <- 1:nlevels(y)
Ylearn <- y[learnind]
## Ylearn <- as.numeric(Ylearn)-1
Xlearn <- X[learnind,]
set.seed(seed)
output.rf <- randomForest(Xlearn, y=Ylearn,...)
Xtest <- X[-learnind,,drop=FALSE]
if(nrow(Xtest) == 0){ Xtest <- Xlearn ; y <- as.numeric(Ylearn)-1 }
else y <- as.numeric(y[-learnind])-1
if(is.null(colnames(Xlearn))) colnames(Xlearn) <- as.character(1:ncol(Xlearn))
colnames(Xtest) <- colnames(Xlearn)
pred.test <- predict(output.rf, newdata=Xtest, type="prob")
yhat <- (0:(nlevels(Ylearn)-1))[apply(pred.test, 1, which.max)]
if(varimp){
	
	modd<-list(NULL)
	if(models==TRUE)
		modd<-list(output.rf)
	
varsel <- as.numeric(importance(output.rf))
new("clvarseloutput", y=y, yhat=yhat, learnind = learnind,
     prob = pred.test, method = "rf", mode=mode, varsel=varsel,model=modd)
}
else{
	modd<-list(NULL)
	if(models==TRUE)
		modd<-list(output.rf)
new("cloutput", y=y, yhat=yhat, learnind = learnind,
     prob = pred.test, method = "rf", mode=mode,model=modd)
}
})

### signature X=matrix, y=factor, f=missing:

setMethod("rfCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, varimp=TRUE, seed=111, models=FALSE,...){
rfCMA(X, y=as.numeric(y)-1, learnind=learnind, varimp=varimp, seed=seed, models=models,...)
})

### signature X=data.frame, f=formula

setMethod("rfCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, varimp=TRUE, seed=111,models=FALSE,...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
rfCMA(as.matrix(X), y=y, learnind=learnind, varimp = varimp, seed = seed, models=models, ...)})


### signature: X=ExpressionSet, y=character.

setMethod("rfCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind, varimp=TRUE, seed=111, models=FALSE, ...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          rfCMA(X=X, y=y, learnind=learnind, varimp=varimp, seed=seed, models=models,...)})