### filename: pknnCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 16.10.2007
#
### Brief description:
#
#  "pknn" stands for probabilisitic nearest neighbour.
#  distances are used as weights and for assigning
#  probabilities. parameter beta regulates the slope of the
#  softmax function.
#
### Further comments and notes:
#
###**************************************************************************###

setGeneric("pknnCMA", function(X, y, f, learnind, beta = 1, k=1, ...)
           standardGeneric("pknnCMA"))

### X=matrix, y=numeric, f=missing

setMethod("pknnCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, beta = 1, k=1, ...){
require(class, quietly=TRUE)
nrx <- nrow(X)
ly <- length(y)
if(nrx != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) stop("'learnind' must not be missing \n")
if(length(learnind) > nrx)
stop("length of 'learnind' must be smaller than the number of observations. \n")
if(k > length(learnind))
stop("'k' chosen too large \n")
y <- as.factor(y)
levels(y) <- 1:nlevels(y)
if(nlevels(y) > 2) mode <- "multiclass"
else mode <- "binary"
y <- as.numeric(y)-1
Ylearn <- y[learnind]
Xlearn <- X[learnind,]
Xtest <- X[-learnind,,drop=FALSE]
if(nrow(Xtest) == 0) stop("Test set is required \n")
else{
y <- y[-learnind]
Xarr <- rbind(Xlearn, Xtest)
DD <- as.matrix(dist(Xarr))
nlearn <- nrow(Xlearn)
part <- DD[1:nlearn,(nlearn+1):nrx]
part <- apply(part, 2, function(z){ ind <- order(z)[1:k]; z[-ind] <- 0; z})
freq <- apply(part, 2, function(z) {z[z > 0] <- Ylearn[which(z > 0)]+1; z})
freq <- apply(freq, 2, function(z) { tab <- table(z[z > 0]); nam <- as.numeric(names(tab))
                                     #tab <- tab[order(nam)];
                                     z[z>0] <- 1/tab[match(z[z > 0], nam)]; z})
part <- part * freq
Gt <- t(model.matrix(~as.factor(Ylearn)-1))
decM <- t(Gt %*% part)
prob <- t(apply(decM, 1, function(z) ifelse(z == 0, 0, exp(-beta*z))))
prob <- prob/rowSums(prob)
if(any(!is.finite(prob)))
warning("class probabilities cannot be computed; reduce size of parameter 'beta' \n")
yhat <- apply(prob, 1, which.max)-1
}
new("cloutput", yhat=yhat, y=y, learnind = learnind,
     prob = prob, method = "pknn", mode=mode)
     })

### signature X=matrix, y=factor, f=missing

setMethod("pknnCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, beta = 1, k=1, ...){
pknnCMA(X, y=as.numeric(y)-1, learnind=learnind, beta = beta, k=k, ...)
})

### signature X=data.frame, f=formula

setMethod("pknnCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, beta = 1, k=1, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
pknnCMA(as.matrix(X), y=y, learnind=learnind, beta = beta, k=k, ...)})


### signature: X=ExpressionSet, y=character.

setMethod("pknnCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind, beta = 1, k=1, ...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          pknnCMA(X=X, y=y, learnind=learnind, beta = beta, k=k, ...)})