### filename: pls_rfCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski, adapted from A-L Boulesteix
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 22.10.2007
#
### Brief description:
#
#  Returns an object of class cloutput.
#
### Further comments and notes:
#
#   - plot functionally included
#   - s.also pls_ldaCMA.r
#   - s.also pls_lrCMA.r
#   ... - argument: passed to randomForest
#   tuning ?
#
###**************************************************************************###

### generic

setGeneric("pls_rfCMA", function(X, y, f, learnind, comp = 2*nlevels(as.factor(y)), seed = 111, models=FALSE,...)
           standardGeneric("pls_rfCMA"))

### X=matrix, y=numeric, f=missing

setMethod("pls_rfCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, comp = 2*nlevels(as.factor(y)), seed=111, models=FALSE,...){
require(plsgenomics, quietly = TRUE)
require(randomForest, quietly = TRUE)
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
Xlearn <- X[learnind,,drop=FALSE]
Xtest <- X[-learnind,,drop=FALSE]
output.pls <- pls.regression(Xlearn, transformy(Ylearn+1), ncomp=comp)
data.learn <- scale(Xlearn, scale=FALSE, center=output.pls$meanX)%*%output.pls$R
if(nrow(Xtest) == 0){ data.test <- data.learn ; y <- Ylearn }
else{
y <- y[-learnind]
data.test <- scale(Xtest,scale=FALSE,center=output.pls$meanX)%*%output.pls$R
}
if(is.null(colnames(data.learn))) colnames(data.learn) <- as.character(1:ncol(data.learn))
colnames(data.test) <- colnames(data.learn)
set.seed(seed)
output.rf <- randomForest(x=data.learn, y = factor(Ylearn), ...)
pred.test <- predict(output.rf, newdata=data.test, type="prob")
yhat <- (0:(length(unique(Ylearn))-1))[apply(pred.test, 1, which.max)]

modd<-list(NULL)
if(models==TRUE)
	modd<-list(output.pls)

new("cloutput", yhat=yhat, y=y, learnind = learnind,
     prob = pred.test, method = "pls_rf", mode=mode,model=modd)
})

### signature X=matrix, y=factor, f=missing:

setMethod("pls_rfCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, comp = 2*nlevels(as.factor(y)), seed = 111, models=FALSE,...){
pls_rfCMA(X, y = as.numeric(y)-1, learnind = learnind, comp = comp, seed = seed, models=models,...)
})

### signature X=data.frame, f=formula

setMethod("pls_rfCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, comp = 2*nlevels(as.factor(y)), seed = 111,models=FALSE, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
pls_rfCMA(as.matrix(X), y=y, learnind=learnind, comp = comp, seed = seed, models=models,...)})


### signature: X=ExpressionSet, y=character.

setMethod("pls_rfCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind, comp = 2*nlevels(as.factor(y)), seed = 111,models=FALSE, ...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          pls_rfCMA(X=X, y=y, learnind=learnind, comp = comp, seed = seed,models=models, ...)})