### filename: svmCMA.r
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
#   The default parameter settings are different from those
#   in e1071 due to the particular setting (p >> n paradigm).
#   Default kernel is the standard scalar product (linear kernel)
#   In addition, the nu parametrization is used.
#
###**************************************************************************###

setGeneric("svmCMA", function(X, y, f, learnind,probability,models=FALSE, ...)
           standardGeneric("svmCMA"))

### X=matrix, y=numeric, f=missing 

setMethod("svmCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, probability, models=FALSE,...){
library(e1071, pos = length(search()))
nrx <- nrow(X)
ly <- length(y)
if(missing(probability)) probability<-FALSE
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
Xlearn <- data.frame(X[learnind,])
dotsCall <- substitute(list(...))
ll <- eval(dotsCall)
if(!hasArg(kernel))ll$kernel <- "linear"
ll$type <- "C-classification"
if(!hasArg(cost)) ll$cost <- 100
ll$x <- Xlearn
ll$y <- factor(Ylearn)
ll$probability <- probability
output.svm <- do.call("svm", args = ll)
Xtest <- X[-learnind,,drop=FALSE]
if(nrow(Xtest) == 0){ Xtest <- Xlearn ; y <- Ylearn }
else y <- y[-learnind]

if(probability==TRUE){ 
pred.test <- predict(object=output.svm, newdata = Xtest, probability = TRUE)
prob <- attr(pred.test, "probabilities")[,as.character(0:(length(unique(Ylearn))-1))]
if(is.vector(prob)) prob <- t(prob)

modd<-list(NULL)
if(models==TRUE)
	modd<-list(output.svm)

ret.obj<-new("cloutput", yhat=as.numeric(pred.test)-1, y=y, learnind = learnind,
     prob = prob, method = "svm", mode=mode,model=modd)

}

if(probability==FALSE){
pred.test <- predict(object=output.svm,newdata=Xtest,probability=F)

modd<-list(NULL)
if(models==TRUE)
	modd<-list(output.svm)

ret.obj<-new("cloutput", yhat=as.numeric(pred.test)-1, y=y, learnind = learnind,
     prob = matrix(data = NA, nrow = length(learnind)), method = "svm", mode=mode,model=modd)


}

ret.obj


})

#### signature X=matrix, y=numeric, f=missing

setMethod("svmCMA", signature(X="matrix", y="factor", f="missing"),#!!?
          function(X, y, learnind, probability, models=FALSE,...){
svmCMA(X, y=as.numeric(y)-1, learnind=learnind,probability=probability,models=models,...)
})

### signature X=data.frame, f=formula

setMethod("svmCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind,probability, models=FALSE,...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
svmCMA(as.matrix(X), y=y, learnind=learnind,probability=probability,models=models,...)})


### signature: X=ExpressionSet, y=character.

setMethod("svmCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind,probability,models=FALSE,...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          svmCMA(X=X, y=y, learnind=learnind,probability=probability, models=models,...)})
