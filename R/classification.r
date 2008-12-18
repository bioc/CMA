### Author: M. Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 9.10.2007
#
### Brief description:
#
#   Returns a list of objects of class "cloutput",
#   for various learningsets.
#   "General" function of the package.
#
### Further comments and notes:
#   s. evaluation.r
#   s. GeneSelection.r
#   s. Classifcation.r
#
###**************************************************************************###

### generic

setGeneric("classification", function(X, y, f, learningsets,
            genesel, genesellist = list(), nbgene, classifier,
            tuneres, tuninglist = list(), trace =TRUE, ...) standardGeneric("classification"))

### X=matrix, y=numeric, f=missing

setMethod("classification", signature(X = "matrix", y = "numeric", f = "missing"),
          function(X, y, f, learningsets, genesel, genesellist = list(),
                   nbgene, classifier, tuneres, tuninglist = list(), trace = TRUE, ...){
dotsCall <- substitute(list(...))
ll <- eval(dotsCall)
if(missing(classifier)) stop("argument 'classifier' is missing \n")

if(missing(learningsets)){
  warning("Argument 'learningsets' is missing; set to a row vector with entries '1:nrow(X)' \n")
  learnmatrix <- matrix(1:nrow(X), ncol=nrow(X))
  }
else{
learnmatrix <- learningsets@learnmatrix
if(ncol(learnmatrix) > nrow(X))
  stop("'learningsets' do not match the input data \n")
}

if(missing(genesel)){
 if(!missing(genesellist) && length(genesellist) != 0){
  genesellist$X <- X
  genesellist$y <- y
  if(!missing(learningsets)) genesellist$learningsets <- learningsets
  genesel <- do.call(GeneSelection, args=genesellist)
  }
 }

else{ if(class(genesel) != "genesel") stop("'genesel' must be of class 'genesel' \n")
      ngenes <- ncol(genesel@rankings[[1]])
      nitergenesel <- nrow(genesel@rankings[[1]])
      if(ngenes != ncol(X)) stop("object 'genesel' does not match the input data \n")
      if(nitergenesel != nrow(learnmatrix))
      stop("object 'genesel' does not match 'learningsets' \n")
    }

if(!missing(nbgene)){
 if(nbgene > ncol(X)) stop("'nbgene' greater than the number all genes \n")}
else nbgene <- ncol(X)

if(missing(tuneres)){
if(!missing(tuninglist) && length(tuninglist) != 0){
 tuninglist$X <- X
 tuninglist$y <- y
 tuninglist$classifier <- classifier
 if(!missing(learningsets)) tuninglist$learningsets <- learningsets
 if(!missing(genesel)){ tuninglist$genesel <- genesel ; tuninglist$nbgene <- nbgene }
 if(!is.list(tuninglist$grid)) stop("Invalid specification of 'tuninglist'. Grid must itself be a list \n")
 tuneres <- do.call(tune, args=c(tuninglist, ll))
 }
}

if(!missing(tuneres)){
 if(length(tuneres@tuneres) != nrow(learnmatrix))
 stop("object 'tuneres' does not match 'learningsets' \n")
 if(!grep(tuneres@method, match.fun(classifier)@generic[1], ignore.case =TRUE))
 stop("object 'tuneres' does not match the chosen classifier. \n")
 besthyperpar <- best(tuneres)
 }

cloutlist <- vector(mode="list", length=nrow(learnmatrix))


 if(missing(genesel)){
 if(missing(tuneres)){
 for(i in 1:nrow(learnmatrix)){
 if(trace) cat("iteration", i, "\n")
 cloutlist[[i]] <- do.call(classifier, args=c(list(X=X, y=y, learnind = learnmatrix[i,]), ll))
  }
 }
 else{
 arglist <- expression(c(ll, list(X=X, y=y, learnind=learnmatrix[i,]), besthyperpar[[i]]))
 funtocall <- as.character(substitute(classifier))
 for(i in 1:nrow(learnmatrix)){
 if(trace) cat("iteration", i, "\n")
 cloutlist[[i]] <- do.call(funtocall, args=eval(arglist))
   }
  }
 }
 else{
 ranks <- genesel@rankings
 imps <- genesel@importance
 if(missing(tuneres)){
 arglist <- expression(c(ll, list(X=Xi, y=y, learnind=learnmatrix[i,])))
 funtocall <- as.character(substitute(classifier))
 }
 else{
   arglist <- expression(c(ll, list(X=Xi, y=y, learnind=learnmatrix[i,]), besthyperpar[[i]]))
   funtocall <- as.character(substitute(classifier))
 }
 if(is.element(genesel@method, c("lasso", "elasticnet", "boosting"))){
  if(length(ranks) > 1){
  for(i in 1:nrow(learnmatrix)){
  if(trace) cat("iteration", i, "\n")
  seli <- c()
   for(j in 1:length(ranks)){
    rankj <- ranks[[j]][i,]
    impj <- imps[[j]][i,]
    impj <- impj[impj > 0]
    nbgene <- min(length(impj), nbgene)
    seli <- c(seli, rankj[1:nbgene])
   }
   Xi <- X[,seli,drop=FALSE]
   cloutlist[[i]] <- do.call(funtocall, args=eval(arglist))
   }
  }
 else{
  ranks <- ranks[[1]]
  imps <- imps[[1]]
  for(i in 1:nrow(learnmatrix)){
  if(trace) cat("iteration", i, "\n")
  impi <- imps[i,]
  impi <- impi[impi > 0]
  nbgene <- min(length(impi), nbgene)
  seli <- ranks[i,1:nbgene]
  Xi <- X[,seli,drop=FALSE]
  cloutlist[[i]] <- do.call(funtocall, args=eval(arglist))
  }
 }
}
else{
  if(length(ranks) > 1){
  for(i in 1:nrow(learnmatrix)){
  if(trace) cat("iteration", i, "\n")
  seli <- c()
   for(j in 1:length(ranks)){
    rankj <- ranks[[j]][i,]
    seli <- c(seli,  rankj[1:nbgene])
   }
   seli <- unique(seli)
   Xi <- X[,seli,drop=FALSE]
   cloutlist[[i]] <- do.call(funtocall, args=eval(arglist))
  }
 }
 else{
  ranks <- ranks[[1]]
  imps <- imps[[1]]
  for(i in 1:nrow(learnmatrix)){
  if(trace) cat("iteration", i, "\n")
  seli <- ranks[i,1:nbgene]
  Xi <- X[,seli,drop=FALSE]
  cloutlist[[i]] <- do.call(funtocall, args=eval(arglist))
   }
  }
 }
}

return(cloutlist)

})

### X=matrix, y=factor, f=missing

setMethod("classification", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learningsets, genesel, genesellist = list(), nbgene,
                  classifier, tuneres, tuninglist = list(), trace =TRUE, ...){
classification(X, y=as.numeric(y)-1, learningsets=learningsets,
               genesel = genesel, genesellist = genesellist, nbgene = nbgene,
               classifier = classifier, tuneres = tuneres, tuninglist = tuninglist,
               trace = trace, ...)
})

### X=matrix, y=missing, f=formula

setMethod("classification", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learningsets, genesel, genesellist = list(), nbgene,
                  classifier, tuneres, tuninglist = list(), trace =TRUE, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
classification(as.matrix(X), y=y, learningsets=learningsets,
               genesel = genesel, genesellist = genesellist, nbgene = nbgene,
               classifier = classifier, tuneres = tuneres, tuninglist = tuninglist,
               trace = trace, ...)})

### X=ExpressionSet, y="character", f="missing"

setMethod("classification", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, f, learningsets, genesel, genesellist = list(), nbgene,
                  classifier, tuneres, tuninglist = list(), trace =TRUE, ...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          classification(X=X, y=y, learningsets=learningsets,
               genesel = genesel, genesellist = genesellist, nbgene = nbgene,
               classifier = classifier, tuneres = tuneres, tuninglist = tuninglist,
               trace = trace, ...)})