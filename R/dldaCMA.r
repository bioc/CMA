setGeneric("dldaCMA", function(X, y, f, learnind, ...)
           standardGeneric("dldaCMA"))


### signature X=matrix, y=numeric, f=missing:

setMethod("dldaCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, ...){
nrx <- nrow(X)
ly <- length(y)
if(nrx != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) learnind <- 1:nrx
if(length(learnind) > nrx)
stop("length of 'learnind' must be smaller than the number of observations. \n")
y <- as.factor(y)
K <- nlevels(y)
levels(y) <- 1:K
if(K > 2) mode <- "multiclass"
else mode <- "binary"
y <- as.numeric(y)-1
Ylearn <- y[learnind]
#X <- data.frame(X)
Xlearn <- X[learnind,]
centroids <- matrix(nrow = ncol(X), ncol = K)
variances <- matrix(nrow = ncol(X), ncol = K)
for(k in 1:K){
  indk <- (Ylearn == (k-1))
  centroids[,k] <- muk <- colMeans(Xlearn[indk,,drop=FALSE])
  variances[,k] <- colSums((Xlearn[indk,,drop=FALSE]-muk)^2)
}
variances <- rowSums(variances)/(length(learnind) - K)
priors <- as.numeric(-2*log(table(Ylearn)/sum(table(Ylearn))))
Xtest <- X[-learnind,,drop=FALSE]
if(nrow(Xtest) == 0){ Xtest <- Xlearn ; y <- Ylearn }
else y <- y[-learnind]
#prob <- predict(object=output.dlda, newdata=Xtest, type = "raw")
#browser()
#yhat <- apply(prob, 1, which.max)-1
dist <- t(apply(Xtest, 1, function(z) colSums((z-centroids)^2/variances) + priors))
dist <- dist - rowMeans(dist)
prob <- safeexp(-0.5*dist)
prob <- prob/rowSums(prob)
yhat <- apply(prob, 1, which.max)-1
new("cloutput", yhat=yhat, y=y, learnind = learnind,
     prob = prob, method = "DLDA", mode=mode)
})

### signature X=matrix, y=factor, f=missing:

setMethod("dldaCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, ...){
dldaCMA(X, y=as.numeric(y)-1, learnind=learnind,...)
})

### signature X=data.frame, f=formula

setMethod("dldaCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
dldaCMA(as.matrix(X), y=y, learnind=learnind,...)})


### signature: X=ExpressionSet, y=character.

setMethod("dldaCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind,...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          dldaCMA(X=X, y=y, learnind=learnind, ...)})