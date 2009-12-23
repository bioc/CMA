### filename: ldaCMA.r
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
###**************************************************************************###

### library(Biobase)  --- for class "ExpressionSet"

### generic

setGeneric("ldaCMA", function(X, y, f, learnind, models=FALSE,...)
           standardGeneric("ldaCMA"))


### signature X=matrix, y=numeric, f=missing:

setMethod("ldaCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, models=FALSE,...){
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
output.lda <- lda(Xlearn, grouping=factor(Ylearn),...)
Xtest <- data.frame(X[-learnind,,drop=FALSE])
if(nrow(Xtest) == 0){ Xtest <- Xlearn ; y <- Ylearn }
else y <- y[-learnind]
pred.test <- predict(object=output.lda, newdata=data.frame(Xtest))

modd<-list(NULL)
if(models==TRUE)
	modd<-list(output.lda)

new("cloutput", yhat=as.numeric(pred.test$class)-1, y=y, learnind = learnind,
     prob = pred.test$posterior, method = "LDA", mode=mode,model=modd)
})

### signature X=matrix, y=factor, f=missing:

setMethod("ldaCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, models=FALSE,...){
ldaCMA(X, y=as.numeric(y)-1, learnind=learnind,models=models,...)
})

### signature X=data.frame, f=formula

setMethod("ldaCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, models=FALSE,...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
ldaCMA(as.matrix(X), y=y, learnind=learnind,models=models,...)})


### signature: X=ExpressionSet, y=character.

setMethod("ldaCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind,models=FALSE,...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          ldaCMA(X=X, y=y, learnind=learnind, models=models,...)})