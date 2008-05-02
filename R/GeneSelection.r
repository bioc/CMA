###
### filename: GeneSelection.r
### Title: Various gene selection methods.
###
### Author: M. Slawski, adapted from A-L Boulesteix
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 26.9.2007
#
### Brief description:
#
#  Returns an object of class 'genesel'.
#
### Further comments and notes:
#
###**************************************************************************###

### generic

setGeneric("GeneSelection", function(X, y, f, learningsets,
            method=c("t.test", "welch.test", "wilcox.test", "f.test", "kruskal.test",
                    "limma", "rfe", "rf", "lasso", "elasticnet", "boosting", "golub"), scheme,
                    trace = TRUE, ...)
           standardGeneric("GeneSelection"))

### X=matrix, y=numeric, f=missing

setMethod("GeneSelection", signature(X="matrix", y="numeric", f="missing"),
        function(X, y, f, learningsets, method=c("t.test", "welch.test", "wilcox.test", "f.test", "kruskal.test",
                    "limma", "rfe", "rf", "lasso", "elasticnet", "boosting", "golub"), scheme, trace = TRUE, ...)
{
nrx <- nrow(X)
ly <- length(y)
if(nrx != length(y))
stop("Number of rows of 'X' must agree with length of y \n")
y <- as.factor(y)
n <- length(y)
maxlvl <- nlevels(y)
tempcall <- as.character(match.call())
if(tempcall[4] == "temp") maxlvl <- 2
#levels(y) <- 1:maxlvl
y <- as.numeric(y)-1
if(missing(learningsets)){
 warning("Argument 'learningsets' is missing; set to a row vector with entries '1:nrow(X)' \n")
 learnmatrix <- matrix(1:nrx, ncol=nrx)
 }

else{
  learnmatrix <- learningsets@learnmatrix
  if(ncol(learnmatrix) > nrx)
  stop("'learningsets' do not match the input data \n")
  }

method <- match.arg(method)
if(!is.element(method, eval(formals(GeneSelection)$method)))
stop("Invalid 'method' specified \n")


if(!missing(scheme)){
if(!is.element(scheme,  c("pairwise", "one-vs-all", "multiclass")))
stop("Invalid 'scheme' specified. Must be one of 'pairwise' or 'one-vs-all' \n")
}
if( maxlvl == 2) scheme <- "pairwise"
else{
      if(missing(scheme) & is.element(method, c("kruskal.test", "f.test", "rf", "boosting", "limma")))
      scheme <- "multiclass"
      if(missing(scheme) & !is.element(method, c("kruskal.test", "f.test", "boosting", "limma"))){
        warning("y has more than two levels, but 'scheme' is not specified and
                  a multiclass method is not used; set to 'one-vs-all' \n")
        scheme <- "one-vs-all"
      }
      if(!missing(scheme) && (scheme == "multiclass" & !is.element(method, c("kruskal.test", "f.test", "boosting", "limma")))){
      warning("scheme is 'multiclass' although a multiclass method is not used;
                  set to 'one-vs-all' \n")
       scheme <- "one-vs-all"
      }

}

niter <- nrow(learnmatrix)
p <- ncol(X)
outrankings <- outimportance <- matrix(nrow=niter, ncol=p)

if( maxlvl == 2  | scheme == "multiclass")
  {
  # rownames.varsel<-character(niter)
  rankings <- importance <- matrix(0, niter, p)
  selfun <- switch(method, t.test = ttest,
                           welch.test = welchtest,
                           wilcox.test =  wilcoxtest,
                           f.test = ftest,
                           kruskal.test = kruskaltest,
                           limma = limmatest,
                           rf = rfCMA,
                           rfe = rfe,
                           lasso = LassoCMA,
                           elasticnet = ElasticNetCMA,
                           boosting = compBoostCMA,
                           golub = golubcrit)
  for (i in 1:niter)
   {
   #Xi<-Xi[!is.na(yi),]
   #yi<-yi[!is.na(yi)]
   if(trace) cat("GeneSelection: iteration",i, "\n")
   # rownames.varsel[i]<-paste("iteration",i)
   outp <- selfun(X, y, learnind=learnmatrix[i,], ...)@varsel
   outrankings[i,] <- ord <- order(outp, decreasing = TRUE)
   outimportance[i,] <- outp[ord]
  #colnames(varsel)<-as.vector(sapply(list(1:p),FUN=paste,".gene",sep=""))
  #rownames(varsel)<-rownames.varsel
  }

  colnames(outrankings) <- paste("rank", 1:p, sep="")
  colnames(outimportance) <- paste("gene", ord, sep="")
  rownames(outrankings) <- rownames(outimportance)  <- paste("iter.", 1:niter, sep="")
  rankings <- importance <- list()
  rankings[[1]] <- outrankings
  importance[[1]] <- outimportance
 }

 else{
  if( scheme == "pairwise" )
 {
 rankings <- importance <- vector(mode="list")
 m <- 1
 for (k in 1:max(y))
  {
  for (j in 0:(k-1))
   {
   temp <- t(apply(learnmatrix, 1, function(z){
   tempindc <- which(!is.element(y[z], c(k,j)))
   z[tempindc] <- 0; z}))
   if(any(rowSums(temp) == 0))
   stop("Scheme 'pairwise' cannot be performed; not each
         learning set contains members of all classes \n")
   temp <- new("learningsets", learnmatrix = temp)
   outp <- GeneSelection(X=X, y=y, learningsets = temp, method=method,...)
   rankings[[paste(k, "vs.", j)]] <- outp@rankings[[1]]
   importance[[paste(k, "vs.", j)]] <- outp@importance[[1]]
   m<-m+1
   }
  }
 }

 if (scheme=="one-vs-all")
 {
 rankings <- importance <- vector(mode="list")
 for (k in 0:max(y))
  {
  check <- apply(learnmatrix, 1, function(z) sum(y[z]==k))
  if(any(check < 1))
  stop("Scheme 'one-vs-all' cannot be performed;
  not each learning set contains members of all classes \n")
  outp <- GeneSelection(X=X, y=as.numeric(y==k), learningsets=learningsets, method=method,...)
  rankings[[paste(k, "vs. rest")]] <- outp@rankings[[1]]
  importance[[paste(k, "vs. rest")]] <- outp@importance[[1]]
  }
 }
}
new("genesel", rankings=rankings, importance=importance, method=method,
    scheme=scheme)

})

### X=matrix, y=factor, f=missing

setMethod("GeneSelection", signature(X="matrix", y="factor", f="missing"),
        function(X, y, f, learningsets, method=c("t.test", "welch.test", "wilcox.test", "f.test", "kruskal.test",
                    "limma", "rfe", "rf", "lasso", "elasticnet", "boosting", "golub"), scheme, trace = TRUE, ...)
         GeneSelection(X, y=as.numeric(y)-1, learningsets=learningsets,
                       method=method, scheme=scheme, trace=trace, ...))

### X=data.frame, y=missing, f=formula

setMethod("GeneSelection", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learningsets, method=c("t.test", "welch.test", "wilcox.test", "f.test", "kruskal.test",
                    "limma", "rfe", "rf", "lasso", "elasticnet", "boosting", "golub"), scheme, trace = TRUE, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
GeneSelection(X=as.matrix(X), y=y, learningsets=learningsets, method=method,
              scheme = scheme, trace = trace, ...)})

### X=ExpressionSet, y=character, f=missing

setMethod("GeneSelection", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learningsets, method=c("t.test", "welch.test", "wilcox.test", "f.test", "kruskal.test",
                   "limma", "rfe", "rf", "lasso", "elasticnet", "boosting", "golub"), scheme,
                   trace = trace, ...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          GeneSelection(X=X, y=y, learningsets=learningsets, method = method,
          scheme = scheme, trace = trace, ...)})