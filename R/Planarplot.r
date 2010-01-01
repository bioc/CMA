### filename: Planarplot.r
### Title: Visualization function.
###
### Author: M. Slawski.
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 3.10.2007
#
### Brief description:
#
#  Plots separation of two or more classes in the plane.
#
### Further comments and notes:
#
#   - Contouring for the multiclass would be useful.
#   - image plot would perhaps be an alternative.
#   - up to now classification is based only on predind.
#     Classification on more + subspace plot would also
#     be interesting.
#
###**************************************************************************###

### generic

setGeneric("Planarplot", function(X, y, f, learnind, predind, classifier, gridsize=100, ...)
           standardGeneric("Planarplot"))

### X=matrix, y=numeric, f=missing

setMethod("Planarplot", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, learnind, predind, classifier, gridsize=100,...){
nrx <- nrow(X)
ly <- length(y)
if(nrx != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(length(predind) != 2)
stop("Exactly two predictors have to be specified \n")
if(any(!is.element(predind, 1:ncol(X))))
stop("Invalid values for 'predind' specified \n")
if(missing(learnind)) learnind <- 1:nrx
if(length(learnind) > nrx)
stop("length of 'learnind' must be smaller than the number of observations. \n")
gridsize <- floor(gridsize)
if(gridsize < 30){
warning("gridsize too low; set to 30. \n")
gridsize <- 30
}
ran1 <- range(X[,predind[1]])
ran1[1] <- ifelse(ran1[1] > 0, ran1[1]*0.8, ran1[1]*1.2)
ran1[2] <- ifelse(ran1[2] > 0, ran1[2]*1.2, ran1[2]*0.8)
ran2 <- range(X[,predind[2]])
ran2[1] <- ifelse(ran2[1] > 0, ran2[1]*0.8, ran2[1]*1.2)
ran2[2] <- ifelse(ran2[2] > 0, ran2[2]*1.2, ran2[2]*0.8)
gridx1 <- seq(from=ran1[1], to=ran1[2], length=gridsize)
gridx2 <- seq(from=ran2[1], to=ran2[2], length=gridsize)
augX <- as.matrix(expand.grid(gridx1, gridx2))
augX <- augX[order(augX[,1], augX[,2]),]
Xpred <- X[,predind]
colnames(augX) <- colnames(Xpred)
XX <- rbind(Xpred, augX)
yy <- c(y, rep(NA, gridsize^2))
outp <- classifier(XX, yy, learnind=learnind, ...)
ll <- ly+1-length(learnind)
uu <- ly+gridsize^2-length(learnind)
predgrid <- outp@yhat[ll:uu]
probgrid <- outp@prob[ll:uu,,drop=FALSE]
plot(augX[,1], augX[,2], xlab = "predictor 1", ylab = "predictor 2", col=predgrid+2,
     main = paste("method:", outp@method, "      learning set: squares     test set: triangles"), cex=0.5, cex.main=1)
points(X[learnind,predind[1]], X[learnind,predind[2]], col=y[learnind]+2, pch=15)
if(length(y[-learnind]) > 0)
points(X[-learnind,predind[1]], X[-learnind,predind[2]], col=y[-learnind]+2, pch=17)
#lunique <- length(unique(y))
#if(lunique == 2)
#contour(gridx1, gridx2, matrix(probgrid[,1], nrow=gridsize, byrow=TRUE), levels=0.5,
#         lwd=2, add=TRUE, drawlabels=FALSE)
#else{
  #for(j in 1:lunique){
  #  intersection <- apply(probgrid, 1, function(z) prod(z[j]-z[-j])+  sum(z[j]-z[-j]))
  #  contour(gridx1, gridx2, matrix(intersection, nrow=gridsize, byrow=TRUE),
  #         lwd=2, levels=0, add=TRUE, drawlabels=FALSE)
  #}
  #for(j in 1:(lunique-1)){
  # for(k in (j+1):lunique){
  #  cat("j:", j, "k:", k, "\n")
  #  todraw <- apply(probgrid, 1, function(z) z[j]-z[k] + sum(z[j] < z[-c(j,k)])+
  #                  sum(z[k] < z[-c(j,k)]))
  #  contour(gridx1, gridx2, matrix(todraw, nrow=gridsize, byrow=TRUE),
  #         lwd=2, levels=0, drawlabels=FALSE, add=TRUE)
  #  }
  # }
#}
})

### X=matrix, y=factor, f=missing

setMethod("Planarplot", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, predind, classifier, gridsize=100, ...){
Planarplot(X, y=as.numeric(y)-1, learnind=learnind, predind=predind, classifier=classifier,
      gridsize=gridsize,...)
})

### signature X=data.frame, f=formula

setMethod("Planarplot", signature(X="data.frame", y="missing", f="formula"),
          function(X, y,f, learnind, predind, classifier, gridsize=100, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
Planarplot(as.matrix(X), y=y, learnind=learnind, predind=predind, classifier=classifier,
      gridsize=gridsize,...)  })

### signature X=ExpressionSet, y=character.

setMethod("Planarplot", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind, predind, classifier, gridsize=100, ...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          Planarplot(X, y=as.numeric(y)-1, learnind=learnind, predind=predind, classifier=classifier,
          gridsize=gridsize,...)})
