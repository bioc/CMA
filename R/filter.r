###
### filename: filter.r
### Title: Gene selection (filter) methods.
###
### Author: M. Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 28.9.2007
#
### Brief description:
#
#  Returns an object of class 'genesel'.
#
### Further comments and notes:
#
#   Are usually not called directly by the user, but via
#   'GeneSelection'.
#
###**************************************************************************###

# [1] ordinary t-test

ttest <- function(X, y, learnind,...){
yprime <- y[learnind]
Xprime <- X[learnind,,drop=FALSE]
yprime <- as.factor(yprime)
levels(yprime) <- 1:nlevels(yprime)
yprime <- as.numeric(yprime)-1
cl0  <-  yprime == 0
xx1 <- Xprime[cl0,,drop=FALSE]
xx2 <- Xprime[!cl0,,drop=FALSE]
df1 <- sum(cl0)-1
df2 <- sum(!cl0)-1
m1 <- colMeans(xx1)
m2 <- colMeans(xx2)
mu1 <- matrix(rep(m1, df1+1), nrow=df1+1, byrow=TRUE)
mu2 <- matrix(rep(m2, df2+1), nrow=df2+1, byrow=TRUE)
ss1 <- colSums((xx1-mu1)^2)
ss2 <- colSums((xx2-mu2)^2)
s <- sqrt((ss1+ss2)/(df1+df2))
se <- s*sqrt(1/(df1+1)+ 1/(df2+1))
stat <- (m1-m2)/se
#if(pvalues) pvals <- 1-pf(ttest^2, 1, df1+df2)
#else pvals <- rep(NA, nrow(x))
new("varseloutput", varsel=abs(stat))
}

### [2] Welch t-test

welchtest <- function(X, y, learnind,...){
yprime <- y[learnind]
Xprime <- X[learnind,,drop=FALSE]
yprime <- as.factor(yprime)
levels(yprime) <- 1:nlevels(yprime)
yprime <- as.numeric(yprime)-1
cl0  <-  yprime == 0
xx1 <- Xprime[cl0,,drop=FALSE]
xx2 <- Xprime[!cl0,,drop=FALSE]
m1 <- colMeans(xx1)
m2 <- colMeans(xx2)
df1 <- sum(cl0)-1
df2 <- sum(!cl0)-1
mu1 <- matrix(rep(m1, df1+1), nrow=df1+1, byrow=TRUE)
mu2 <- matrix(rep(m2, df2+1), nrow=df2+1, byrow=TRUE)
ss1 <- colSums((xx1-mu1)^2)
ss2 <- colSums((xx2-mu2)^2)
se1 <- ss1/(df1*(df1+1))
se2 <- ss2/(df2*(df2+1))
se12 <-  sqrt(se1 + se2)
stat <- (m1-m2)/se12
new("varseloutput", varsel=abs(stat))
}

### [3] Wilcoxon-Test

wilcoxtest <- function(X, y, learnind,...){

yprime <- y[learnind]
Xprime <- X[learnind,,drop=FALSE]
yprime <- as.factor(yprime)
levels(yprime) <- 1:nlevels(yprime)
yprime <- as.numeric(yprime)-1
cl0  <-  yprime == 0
n0 <- sum(cl0)
n1 <- sum(!cl0)
Rx <- apply(Xprime, 2, rank)
stat <- colSums(Rx[cl0, ,drop=FALSE])-sum(1:sum(cl0))
new("varseloutput", varsel=stat)
}

### [4] F-Test ~ 5.sec 4000 genes, 40 samples

ftest <- function(X, y, learnind,...){
yprime <- y[learnind]
Xprime <- X[learnind,,drop=FALSE]
yprime <- as.factor(yprime)
ll <- length(yprime)
dfnum <- nlevels(yprime)-1
dfdenom <- ll-dfnum-1
colffun <- function(z){
 ssfull <- sum(tapply(z, yprime, function(w) sum((w-mean(w))^2)))
 sstot <- var(z)*(ll-1)
 ((sstot-ssfull)/dfnum )/(ssfull/dfdenom)
}
stat <- apply(Xprime, 2, function(z) colffun(z))
new("varseloutput", varsel=stat)
}

### [5] http://en.wikipedia.org/wiki/Kruskal-Wallis_test ~5 sec.

kruskaltest <- function(X, y, learnind,...){
yprime <- y[learnind]
Xprime <- X[learnind,,drop=FALSE]
yprime <- as.factor(yprime)
ll <- length(yprime)
const1 <- 12/(ll*(ll+1))
const2 <- (ll+1)/2
Rx <- apply(Xprime, 2, rank)
colfun <- function(z) sum(tapply(z, yprime, function(w) length(w)*(mean(w)-const2)^2))
stat <- apply(Rx, 2, function(z) colfun(z))
stat <- const1 * stat
new("varseloutput", varsel=stat)
}

### [6] limmatest

limmatest <- function(X, y, learnind,...){
require(limma, quietly=TRUE)
yprime <- y[learnind]
Xprime <- X[learnind,,drop=FALSE]
yprime <- as.factor(yprime)
des <- model.matrix(~yprime)
limo <- lmFit(t(Xprime), des)
cont<-matrix(0,nrow=ncol(des)-1,ncol=ncol(des))
cont[,-1]<-diag(ncol(des)-1)
conp<-contrasts.fit(limo,contrasts=t(cont))
outp <- eBayes(conp, ...)
stat <-  classifyTestsF(object = outp, fstat.only = TRUE)
attributes(stat) <- NULL

new("varseloutput", varsel=stat)
}

### [7] golub criterion

golubcrit <-function(X,y, learnind,...)
{
yprime <- y[learnind]
Xprime <- X[learnind,,drop=FALSE]
cl0  <-  yprime == 0
xx1 <- Xprime[cl0,,drop=FALSE]
xx2 <- Xprime[!cl0,,drop=FALSE]
m1 <- colMeans(xx1)
m2 <- colMeans(xx2)
df1 <- sum(cl0)-1
df2 <- sum(!cl0)-1
sd1 <- sqrt(colSums((xx1-m1)^2)/df1)
sd2 <- sqrt(colSums((xx2-m2)^2)/df2)
stat <- (m1 - m2)/(sd1 + sd2)
new("varseloutput", varsel=abs(stat))
}

### [8] Recursive feature elimination

rfe <- function(X,y,learnind,...){
require(e1071, quietly=TRUE)
dots <- eval(substitute(list(...)))
if(!hasArg(kernel)) dots$kernel <- "linear"
dots$type <- "C-classification"
if(!hasArg(cost)) dots$cost <- 100
dots$y <- y[learnind]
dots$x <- X[learnind,,drop=FALSE]
model <- do.call(svm, args=dots)
stat <- as.numeric((t(model$coefs) %*% dots$x[model$index,,drop=FALSE])^2)
new("varseloutput", varsel=stat)
}

### [9] shrinked correlation adjusted t ("shrinkcat")
shrinkcat<-function(X,y,learnind,...){
 ##library
 require('st', quietly=TRUE)
 ## arguments for shrinkcat.stat call
 if(nlevels(y)>2) stop("The class vector must have 2 labels")	
 dots<-eval(substitute(list(...)))
 dots$X=X[learnind,,drop=FALSE]
 dots$L=y[learnind]
 dots$verbose=FALSE
 stat<-do.call('shrinkcat.stat', args=dots)
 new("varseloutput", varsel=abs(stat))
}









