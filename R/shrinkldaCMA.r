### filename: shrinkldaCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 18.9.2007
#
### Brief description:
#
#
#  Returns an object of class cloutput.
#
### Further comments and notes:
#
#   Requires package 'corpcor'.
#
###**************************************************************************###

setGeneric("shrinkldaCMA", function(X, y, f, learnind, ...)
           standardGeneric("shrinkldaCMA"))


### signature X=matrix, y=numeric, f=missing:

setMethod("shrinkldaCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, ...){
require(corpcor, quietly=TRUE)
nrx <- nrow(X)
ly <- length(y)
if(nrow(X) != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) learnind <- 1:nrx
y <- as.factor(y)
levels(y) <- 1:nlevels(y)
if(nlevels(y) > 2) mode <- "multiclass"
else mode <- "binary"
y <- as.numeric(y)-1
Ylearn <- y[learnind]
Xlearn <- X[learnind,,drop=FALSE]
taby <- table(Ylearn)
priors <- as.numeric(taby/sum(taby))
svdXlearn <- svd(Xlearn)
svalue <- svdXlearn$d
svaluePos <- seq(svalue)[svalue > 0]
svalue <- svalue[svaluePos]
R <- svdXlearn$u[, svaluePos] %*% diag(svalue)
G <- model.matrix(~factor(Ylearn) - 1)
centroids <-  scale(t(R) %*% G, FALSE, taby)
Sigmastar <- cov.shrink(R, verbose=FALSE, ...)
Sigmastarinv <- solve(Sigmastar)
coefs <- svdXlearn$v %*% Sigmastarinv %*% centroids
QF <- apply(centroids, 2, function(z) -0.5*mahalanobis(z, center=rep(0, length(centroids)),
                                               cov=Sigmastarinv, inverted=TRUE))
Xtest <- X[-learnind,,drop=FALSE]
if(nrow(Xtest) == 0){ Xtest <- Xlearn ; y <- Ylearn }
else y <- y[-learnind]
Dis <-  Xtest %*% coefs + QF + log(priors)
classes <- as.numeric(names(taby)[apply(Dis, 1, which.max)])
probs <- exp(Dis)
probs <- probs/rowSums(probs)
colnames(probs) <- names(taby)

new("cloutput", y=y, yhat=classes, learnind = learnind,
    prob = probs, method = "shrinkLDA", mode=mode)
})

### signature X=matrix, y=factor, f=missing

setMethod("shrinkldaCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, ...){
shrinkldaCMA(X, y=as.numeric(y)-1, learnind=learnind,...)
})

### signature X=data.frame, f=formula

setMethod("shrinkldaCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
shrinkldaCMA(as.matrix(X), y=y, learnind=learnind,...)})


### signature: X=ExpressionSet, y=character.

setMethod("shrinkldaCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind,...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          shrinkldaCMA(X=X, y=y, learnind=learnind, ...)})