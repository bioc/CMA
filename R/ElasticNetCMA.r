### filename: ElasticNetCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 25.9.2007
#
### Brief description:
#
#  Returns an object of class clvarseloutput.
#
### Further comments and notes:
#
#   Wrapper to the glmpath package.
#   s. ElasticNetCMA.r
#
#   involves a second hyperparameter.
#
###**************************************************************************###

setGeneric("ElasticNetCMA", function(X, y, f, learnind, norm.fraction = 0.1, lambda2=1e-3, ...)
           standardGeneric("ElasticNetCMA"))


setMethod("ElasticNetCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, norm.fraction=0.1, lambda2=1e-3, ...){
require(glmpath, quietly=TRUE)
nrx <- nrow(X)
ly <- length(y)
if(nrow(X) != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) learnind <- 1:nrx
y <- as.factor(y)
levels(y) <- 1:nlevels(y)
y <- as.numeric(y)-1
Ylearn <- y[learnind]
if(nlevels(Ylearn) > 2) stop("'ElasticNetCMA' only possible for binary classification \n")
mode <- "binary"
Xlearn <- X[learnind,]
ll <- eval(substitute(list(...)))
ll$family <- binomial
ll$lambda2 <- lambda2
ll$x <- Xlearn
ll$y <- Ylearn
output <- do.call(glmpath, args=ll)
norm.fraction <- norm.fraction[1]
varsel <- predict(output, s=norm.fraction, mode="norm.fraction", type="coefficients")
varsel <- abs(as.numeric(varsel[-1]))
Xtest <- X[-learnind,,drop=FALSE]
if(nrow(Xtest) == 0) { Xtest <- Xlearn ; y <- Ylearn }
else y <- y[-learnind]
prob <- as.numeric(predict(output, newx=Xtest, s=norm.fraction, mode="norm.fraction",
                   type="response"))
yhat <- as.numeric(prob > 0.5)
prob <- cbind(1-prob, prob)

new("clvarseloutput", y=y, yhat=yhat, learnind = learnind,
     prob = prob, method = "ElasticNet", mode=mode, varsel=varsel)
})

### signature X=matrix, y=factor, f=missing:

setMethod("ElasticNetCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, norm.fraction=0.1, lambda2=1e-3, ...){
ElasticNetCMA(X, y=as.numeric(y)-1, learnind=learnind,
              norm.fraction = norm.fraction, lambda2 = lambda2, ...)
})

### signature X=data.frame, f=formula

setMethod("ElasticNetCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, norm.fraction=0.1, lambda2=1e-3, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
ElasticNetCMA(as.matrix(X), y=y, learnind=learnind, norm.fraction = norm.fraction,
              lambda2 = lambda2, ...)})


### signature: X=ExpressionSet, y=character.

setMethod("ElasticNetCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind, norm.fraction=0.1, lambda2=1e-3, ...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          ElasticNetCMA(X=X, y=y, learnind=learnind, norm.fraction = norm.fraction,
              lambda2 = lambda2, ...)})