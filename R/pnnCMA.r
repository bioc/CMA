### filename: pnnCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 18.9.2007
#
### Brief description:
#  In fact a Nadaraya-Watson kernel estimator (gaussian kernel), re-labeled
#  as a probabilistic neural network by Specht (1990, Neural Networks).
#
#  Returns an object of class cloutput.
#
### Further comments and notes:
#   classification performance is sometimes very bad and very unstable,
#   might profit from variable selection. Data are scaled before
#   training.
#
###**************************************************************************###

### generic

setGeneric("pnnCMA", function(X, y, f, learnind, sigma = 1,models=FALSE)
           standardGeneric("pnnCMA"))

### signature X=matrix, y=numeric, f=missing

setMethod("pnnCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, learnind, sigma=1,models=FALSE){
X <- scale(X)
X <- as.matrix(X)
nrx <- nrow(X)
ly <- length(y)
if(nrx != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) stop("'learnind' must not be missing \n")
if(length(learnind) > nrx)
stop("length of 'learnind' must be smaller than the number of observations. \n")
y <- as.factor(y)
layers <- nlevels(y)
levels(y) <- 1:layers
if(layers > 2) mode <- "multiclass"
else mode <- "binary"
y <- as.numeric(y)-1
Ylearn <- y[learnind]
W <- X[learnind,,drop=FALSE]
V <- X[-learnind,,drop=FALSE]
if(nrow(V) == 0) stop("Test set is required \n")
Wlist <- sapply(0:(layers-1), function(z) subset(W, Ylearn==z))
rbffun <- function(z) dnorm(z, mean = 0, sd = sigma)
res <- apply(V, 1, function(v){
                         summation <- numeric(layers)
                         for(i in 1:layers){
                          summation[i] <- mean(rbffun(rowSums((Wlist[[i]]- v)^2)))
                          #summation[i] <- mean(rbffun(colSums(v%*% t(Wlist[[i]]) - 1)))
                          }
                          return(summation)
                        })
res <- t(res)
if(any(rowSums(res) == 0))
warning("Value for 'sigma' should be changed, all output units  are zero for
      at least one test observation; NaNs for class probabilities occured \n")
yhat <- (0:(layers-1))[apply(res, 1, which.max)]
prob <- res/rowSums(res)

if(models==TRUE)
	modd<-list(NULL)
if(models==FALSE)
	modd<-list(NULL)

new("cloutput", yhat=yhat, y=y[-learnind], learnind = learnind,
     prob = prob, method = "pnn", mode=mode,model=modd)
})

### signature X=matrix, y=factor, f=missing

setMethod("pnnCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, sigma = 1,models=FALSE){
pnnCMA(X, y=as.numeric(y)-1, learnind=learnind, sigma = sigma,models=models)
})

### signature X=data.frame, f=formula

setMethod("pnnCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, sigma = 1,models=FALSE){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
pnnCMA(as.matrix(X), y=y, learnind=learnind, sigma = sigma,models=models)})


### signature: X=ExpressionSet, y=character.

setMethod("pnnCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind, sigma = 1,models=FALSE){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          pnnCMA(X=X, y=y, learnind=learnind, sigma=sigma,models=models)})