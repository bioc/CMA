### filename: LassoCMA.r
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
###**************************************************************************###

setGeneric("LassoCMA", function(X, y, f, learnind, norm.fraction = 0.1, ...)
           standardGeneric("LassoCMA"))


setMethod("LassoCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, norm.fraction=0.1, ...){
require(glmnet, quietly=TRUE)
#require(glmpath) #vgl
nrx <- nrow(X)
ly <- length(y)
if(nrow(X) != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) learnind <- 1:nrx
y <- as.factor(y)
levels(y) <- 1:nlevels(y)
if(nlevels(y) > 2) stop("'LassoCMA' only possible for binary classification \n")
y <- as.numeric(y)-1
Ylearn <- y[learnind]
mode <- "binary"   
type<-'naive'

### levels(Ylearn) <- 1:nlevels(Ylearn)
Xlearn <- X[learnind,]
ll <- eval(substitute(list(...)))
ll$family <- 'binomial'
ll$alpha <- 1
ll$x <- Xlearn
ll$y <- Ylearn

##lambdamin
if(is.null(ll$lambda.min))
  ll$lambda.min <- 1e-03

output <- do.call("glmnet", args=ll)
norm.fraction <- norm.fraction[1]

##determine lambda that corresponds to specified norm.fraction


 nbeta=rbind2(output$a0,output$beta)
nbeta <- as.matrix(nbeta)
b<-t(nbeta)
s <- norm.fraction
k <- nrow(b)
    std.b <- scale(b[ ,-1], FALSE, 1/sd(Xlearn))
 bnorm <- apply(abs(std.b), 1, sum)
   sb <- bnorm / bnorm[k]
sfrac <- (s - sb[1])/(sb[k] - sb[1])
    sb <- (sb - sb[1])/(sb[k] - sb[1])
    usb <- unique(sb)
    useq <- match(usb, sb)
    sb <- sb[useq]
    b <- b[useq, ]
coord <- approx(sb, seq(sb), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    newb <- ((sb[right] - sfrac) * b[left, , drop = FALSE] + 
        (sfrac - sb[left]) * b[right, , drop = FALSE])/(sb[right] - sb[left])
    newb[left == right, ] <- b[left[left == right], ]
betann3<-newb
###########################

varsel <- newb
varsel <- abs(as.numeric(varsel[-1]))
Xtest <- X[-learnind,,drop=FALSE]
if(nrow(Xtest) == 0) { Xtest <- Xlearn ; y <- Ylearn }
else y <- y[-learnind]
prob <- as.numeric(predict(output, newx=Xtest, s=(output$lambda[left]+output$lambda[right])/2,
                   type="response"))
yhat <- as.numeric(prob > 0.5)
prob <- cbind(1-prob, prob)
output.net <- output


new("clvarseloutput", y=y, yhat=yhat, learnind = learnind,
     prob = prob, method = "Lasso", mode=mode, varsel=varsel)
})

### signature X=matrix, y=factor, f=missing:

setMethod("LassoCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, norm.fraction=0.1, ...){
LassoCMA(X, y=as.numeric(y)-1, learnind=learnind, norm.fraction=norm.fraction, ...)
})

### signature X=data.frame, f=formula

setMethod("LassoCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, norm.fraction=0.1, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
LassoCMA(as.matrix(X), y=y, learnind=learnind, norm.fraction=norm.fraction, ...)})


### signature: X=ExpressionSet, y=character.

setMethod("LassoCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind, norm.fraction = 0.1, ...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          LassoCMA(X=X, y=y, learnind=learnind, norm.fraction=norm.fraction, ...)})

      

