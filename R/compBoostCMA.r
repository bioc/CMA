setGeneric("compBoostCMA", function(X, y, f, learnind, loss=c("binomial", "exp", "quadratic"),
                   mstop=100, nu=0.1, models=FALSE,...)
           standardGeneric("compBoostCMA"))


### signature X=matrix, y=numeric, f=missing:

setMethod("compBoostCMA", signature(X="matrix", y="numeric", f="missing"),
          function(X, y, f, learnind, loss=c("binomial", "exp", "quadratic"),
                   mstop=100, nu=0.1, models=FALSE, ...){
nrx <- nrow(X)
ly <- length(y)
if(nrow(X) != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
if(missing(learnind)) learnind <- 1:nrx
y <- as.factor(y)
levels(y) <- 1:nlevels(y)
loss <- match.arg(loss)
if(nlevels(y) > 2){
mode <- "multiclass"
  if(loss != "binomial"){
  warning("For the multiclass case, only the binomial loss is implemented.")
  loss <- "binomial"
  }
}
else mode <- "binary"
y <- as.numeric(y)-1
Ylearn <- y[learnind]
levels(Ylearn) <- 1:nlevels(Ylearn)
Xlearn <- X[learnind,]
if(missing(mstop)) mstop <- 100
if(missing(nu)) nu <- 0.1
if(mode == "binary"){
Ylearn <- 2*(Ylearn)-1
if(!is.element(loss, c("binomial", "exp", "quadratic")))
stop("'loss' must be one of 'binomial', 'exp', 'quadratic'")
ngradient <- switch(loss, binomial = function (y, ff){
                                     exp2yf <- exp(-2 * y * ff)
                                    -(-2 * y * exp2yf)/(log(2) * (1 + exp2yf))},
                          exp= function(y,ff) y * exp(-y * ff),
                          quadratic= function(y,ff){
                           f  <- sign(ff) * pmin(abs(ff), 1)
                           -2 * y + 2 * y * ff })
varsel <- rep(0, ncol(Xlearn))

denom <- colSums(Xlearn^2)
part  <- t(Xlearn)/sqrt(denom)
if (all(is.na(part))) warning("Column-wise inverses cannot be computed \n")
pi0 <- mean(Ylearn > 0)
offset <- 0.5*log(pi0/(1-pi0))
fit <- rep(offset, length(Ylearn))
for (m in 1:mstop){
u <- ngradient(Ylearn, fit)
xselect <- which.max(abs(coef <- part %*% u))
coef <- coef[xselect]/sqrt(denom[xselect])
fit <- fit + (nu * coef) * Xlearn[, xselect]
varsel[xselect] <- varsel[xselect] + coef*nu
}
Xtest <- X[-learnind,,drop=FALSE]
if(nrow(Xtest) == 0){ Xtest <- Xlearn ; y <- (Ylearn + 1)/2}
else y <- y[-learnind]
varselind <- which(varsel != 0)
pred.test <- Xtest[,varselind,drop=FALSE] %*% varsel[varselind] + offset
prob <- plogis(pred.test)
yhat <- as.numeric(pred.test > 0)
prob <- cbind(1-prob, prob)
varsel <- abs(varsel)
}
if(mode == "multiclass"){
 G <- model.matrix(~factor(Ylearn) - 1)
 G <- G[1:nrow(G),1:ncol(G),drop=FALSE]
 colnames(G) <- rownames(G) <- NULL
 ngradient <- function(ff){
  Ef <- exp(ff)
  denom <- rowSums(Ef)
  normEf <- Ef/denom
  U <- G - normEf
  return(U)
 }
 freq <- colSums(G)
 freq <- freq/sum(freq)
 offset <- numeric(ncol(G))
 for(j in seq(along=freq)) offset[j] <- log(freq[j]/(1-freq[j]))
 offset <- offset - mean(offset)
 fit <-  matrix(data=rep(offset, times=length(Ylearn)), ncol=ncol(G), byrow=TRUE)
 varsel <- matrix(nrow=ncol(Xlearn), ncol=ncol(G), data=0)
 varselind <- numeric(mstop)
 sx2 <- sqrt(colSums(Xlearn^2))
 for(m in 1:mstop){
 U <- ngradient(fit)
 B <- crossprod(Xlearn, U)/sx2
 B2 <- rowSums(abs(B))
 xselect <- varselind[m] <- which.max(B2)
 coef <- B[xselect,,drop=FALSE]/sx2[xselect]
 coef <- coef - mean(coef)
 varsel[xselect,] <- varsel[xselect,] + coef*nu
 fit <- fit + nu *  Xlearn[,xselect,drop=FALSE] %*% coef
 }
 Xtest <- X[-learnind,,drop=FALSE]
 if(nrow(Xtest) == 0){ Xtest <- Xlearn  ; y <- Ylearn }
 else y <- y[-learnind]
 varselind <- setdiff(unique(varselind),0)
 pred.test <- Xtest[,varselind,drop=FALSE] %*% varsel[varselind,,drop=FALSE] + offset
 pred.test <- exp(pred.test)
 prob <- pred.test/rowSums(pred.test)
 yhat <- apply(pred.test, 1, which.max)-1
 varsel <- rowSums(abs(varsel))
}

if(models==TRUE)
	modd<-list(NULL)
if(models==FALSE)
	modd<-list(NULL)


new("clvarseloutput", y=y, yhat=yhat, learnind = learnind,
     prob = prob, method = "compBoost", mode=mode, varsel=varsel,model=modd)
})

### signature X=matrix, y=factor, f=missing:

setMethod("compBoostCMA", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learnind, loss=c("binomial", "exp", "quadratic"),
                   mstop=100, nu=0.1, models=FALSE,...){
compBoostCMA(X, y=as.numeric(y)-1, learnind=learnind, loss=loss,
             mstop = mstop, nu=nu, models=models,...)
})

### signature X=data.frame, f=formula

setMethod("compBoostCMA", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learnind, loss=c("binomial", "exp", "quadratic"),
                   mstop=100, nu=0.1, models=FALSE,...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
#f <- as.formula(paste("~", paste(xvar, collapse = "+")))
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
compBoostCMA(as.matrix(X), y=y, learnind=learnind, loss = loss,
              mstop = mstop, nu = nu, models=models, ...)})


### signature: X=ExpressionSet, y=character.

setMethod("compBoostCMA", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learnind, loss=c("binomial", "exp", "quadratic"),
                   mstop=100, nu=0.1,models=FALSE,...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          compBoostCMA(X=X, y=y, learnind=learnind,
                       loss = loss, mstop = mstop, nu = nu,models=models,...)})