### filename: pls_lrCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski
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
#   - s.also pls_rfCMA.r
#
###**************************************************************************###

###

setGeneric("pls_lrCMA", function(X, y, f, learnind, comp = 2, lambda = 1e-4, plot = FALSE)
           standardGeneric("pls_lrCMA"))

### X=matrix, y=numeric, f=missing

setMethod("pls_lrCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, comp = 2, lambda = 1e-4, plot = FALSE){
require(plsgenomics, quietly=TRUE)
nrx <- nrow(X)
ly <- length(y)
if(nrx != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) learnind <- 1:nrx
if(length(learnind) > nrx)
stop("length of 'learnind' must be smaller than the number of observations. \n")
y <- as.factor(y)
levels(y) <- 1:nlevels(y)
if(nlevels(y) > 2) stop("pls_lrCMA is only possible for binary classification \n")
else mode <- "binary"
y <- as.numeric(y)-1
if(plot) yy <- c(y[learnind], y[-learnind])
Ylearn <- y[learnind]
Xlearn <- X[learnind,,drop=FALSE]
Xtest <- X[-learnind,,drop=FALSE]
output.pls <- pls.regression(Xlearn, transformy(Ylearn+1), ncomp=comp)
data.learn <- scale(Xlearn,scale=FALSE,center=output.pls$meanX)%*%output.pls$R
if(nrow(Xtest) == 0){ data.test <- data.learn ; y <- Ylearn }
else{
y <- y[-learnind]
data.test <- scale(Xtest,scale=FALSE,center=output.pls$meanX)%*%output.pls$R
}
design <- cbind(1, data.learn)
output.glm <- penlogitfit(design, y=Ylearn, lambda = lambda)
prob <- plogis(cbind(1, data.test) %*% output.glm)
prob <- cbind(1-prob, prob)
yhat <- apply(prob, 1, which.max)-1

if(plot){
 if(comp > 2) stop("If 'plot=TRUE', number of components have to be 1 or 2. \n")
 else{
    pchs <- ifelse(1:ly %in% learnind, 15, 17)
    if(comp == 1){
      if(nrow(Xtest) == 0){
      plot(data.test[,1], rep(0.5, ly), pch=pchs, col=yy+2, xlab="first pls component",
           main="learning set: square     test set: triangle",
           cex.main=1, ylab="")
      }
      else{
      plot(rbind(data.learn, data.test)[,1], rep(0.5, ly), pch=pchs, col=yy+2, xlab="first pls component",
           main="learning set: square     test set: triangle",
           cex.main=1, ylab="")
      }
    }
    else{
      if(nrow(Xtest) == 0){
      plot(data.learn[,1], data.learn[,2], pch=pchs, col=yy+2, xlab="first pls component",
           main="learning set: square     test set: triangle",
           cex.main=1, ylab="second pls component")
      }
      else{
      dat <- rbind(data.learn, data.test)
      plot(dat[,1], dat[,2], pch=pchs, col=yy+2, xlab="first pls component",
            main="learning set: square     test set: triangle",
            cex.main=1, ylab="second pls component")
      }
    }
   }
  }

new("cloutput", yhat = yhat, y=y, learnind = learnind,
     prob = prob, method = "pls_lr", mode=mode)


})

### signature X=matrix, y=factor, f=missing:

setMethod("pls_lrCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, comp = 2, lambda = 1e-4, plot = FALSE){
pls_lrCMA(X, y = as.numeric(y)-1, learnind = learnind, comp = comp, lambda = lambda, plot = plot)
})

### signature X=data.frame, f=formula

setMethod("pls_lrCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, comp = 2, lambda = 1e-4, plot = FALSE){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
pls_lrCMA(as.matrix(X), y=y, learnind=learnind, comp = comp, lambda = lambda, plot = plot)})


### signature: X=ExpressionSet, y=character.

setMethod("pls_lrCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind, comp = 2, lambda = 1e-4, plot = FALSE){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          pls_lrCMA(X=X, y=y, learnind=learnind, comp = comp, lambda = lambda, plot = plot)})