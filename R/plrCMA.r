### filename: plrCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski, adapted from Ji Zhu
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

setGeneric("plrCMA", function(X, y, f, learnind, lambda = 0.01, scale = TRUE, ...)
           standardGeneric("plrCMA"))

### X=matrix, y=numeric, f=missing

setMethod("plrCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, lambda=0.01, scale =TRUE, ...){
if(scale) X <- scale(X)
nrx <- nrow(X)
ly <- length(y)
if(nrx != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) learnind <- 1:nrx
y <- as.factor(y)
levels(y) <- 1:nlevels(y)
if(nlevels(y) > 2) mode <- "multiclass"
else mode <- "binary"
y <- as.numeric(y)-1
Ylearn <- y[learnind]
Xlearn <- X[learnind,]
gram <- tcrossprod(Xlearn)
if(nrow(X[-learnind,,drop=FALSE]) != 0){
gramtest <- X[-learnind,,drop=FALSE] %*% t(Xlearn)
y <- y[-learnind]
}
else { gramtest <- gram  ; y <- Ylearn }
if(mode == "binary"){
output <- bklr(Ylearn, Ka=gram, Kp=gram, lambda=lambda)
pred.out <- bklr.predict(output$alpha, gramtest)
yhat <- as.numeric(pred.out$fit > 0)
prob <- cbind(1-pred.out$mu, pred.out$mu)
}
if(mode == "multiclass"){
 G <- model.matrix( ~ as.factor(Ylearn)-1)
 k <- length(unique(Ylearn))
 gramk <- array(rep(gram, times=k), dim=c(nrow(gram), nrow(gram), k))
 output <- mklr(G, Ka=gramk, lambda=lambda)
 gramtestk <- array(rep(gramtest, times=k), dim=c(nrow(gramtest), ncol(gramtest), k))
 G <- model.matrix(~ factor(y, levels = levels(as.factor(Ylearn)))-1)
 ### y only used to compute deviance, not for prediction !
 pred.out <- mklr.predict(output, gramtestk, y=G)
 yhat <- pred.out$predict-1
 prob <- pred.out$mu
 }

 #browser()
new("cloutput", y=y, yhat=yhat, learnind = learnind,
     prob = prob, method = "plr", mode=mode)
})

### X=matrix, y=factor, f=missing

setMethod("plrCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, lambda = 0.01, scale =TRUE, ...){
plrCMA(X, y=as.numeric(y)-1, learnind=learnind, lambda = lambda, scale = scale, ...)
})

### signature X=data.frame, f=formula

setMethod("plrCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, lambda = 0.01, scale =TRUE, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
plrCMA(as.matrix(X), y=y, learnind=learnind, lambda = lambda, scale = scale, ...)})


### signature: X=ExpressionSet, y=character.

setMethod("plrCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind, lambda = 0.01, scale = TRUE, ...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          plrCMA(X=X, y=y, learnind=learnind, lambda = lambda, scale = scale, ...)})