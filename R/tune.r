### filename: tune.r
### Title: Function to tune different classication methods.
###
###
### Author: M. Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 9.10.2007
#
### Brief description:
#
#   Returns an object of class "tuningresult",
#   based on results of "evaluation".
#
### Further comments and notes:
#   s. evaluation.r
#   s. GeneSelection.r
#   s. Classifcation.r
#
###**************************************************************************###

### generic

setGeneric("tune", function(X, y, f, learningsets,
            genesel, genesellist = list(), nbgene, classifier, fold = 3, strat = FALSE, grids = list(), trace=TRUE, ...)
           standardGeneric("tune"))

### X=matrix, y=numeric, f=missing

setMethod("tune", signature(X = "matrix", y = "numeric", f = "missing"),
          function(X, y, f, learningsets, genesel, genesellist = list(),
                   nbgene, classifier, fold = 3, strat = FALSE, grids = list(), trace = TRUE, ...){

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
  genesel <- do.call("GeneSelection", args=genesellist)
  }
 }


else{ if(class(genesel) != "genesel") stop("'genesel' must be of class 'genesel' \n")
      ngenes <- ncol(genesel@rankings[[1]])
      nitergenesel <- nrow(genesel@rankings[[1]])
      if(ngenes != ncol(X)) stop("object 'genesel' does not match the input data \n")
      if(nitergenesel != nrow(learnmatrix))
      stop("object 'genesel' does not match 'learningsets' \n")
      }

if(!missing(nbgene))
 if(nbgene > ncol(X)) stop("'nbgene' greater than the number all genes \n")
else nbgene <- ncol(X)

ll <- eval(substitute(list(...)))

classifname <- getMethod(classifier, signature(X="matrix", y="numeric", f="missing"))@generic

if(length(grids) == 0){
 grids <- switch(classifname,#as.character(substitute(classifier)),
                              gbmCMA = list(n.trees = c(50, 100, 200, 500, 1000)),
                              compBoostCMA = list(mstop = c(50, 100, 200, 500, 1000)),
                              LassoCMA = list(norm.fraction = seq(from=0.1, to=0.9, length=9)),
                              ElasticNetCMA = list(norm.fraction = seq(from=0.1, to=0.9, length=5),
                                                   lambda2 = 2^{-(5:1)}),
                              plrCMA = list(lambda = 2^{-4:4}),
                              pls_ldaCMA = list(comp = 1:10),
                              pls_lrCMA = list(comp = 1:10),
                              pls_rfCMA = list(comp = 1:10),
                              rfCMA = list(mtry = ceiling(c(0.1, 0.25, 0.5, 1, 2)*sqrt(ncol(X))),
                                           nodesize = c(1,2,3)),
                              knnCMA = list(k=1:10),
                              pknnCMA = list(k = 1:10),

                              scdaCMA = list(delta = c(0.1, 0.25, 0.5, 1, 2, 5)),
                              pnnCMA = list(sigma = c(2^{-2:2})),
                              nnetCMA = list(size = 1:5, decay = c(0, 2^{-(4:1)})))
 if(classifname == "svmCMA"){
   if(!hasArg(kernel)) ll$kernel <- "linear"
   else ll$kernel <- match.arg(kernel)
    grids <- switch(ll$kernel, linear = list(cost = c(0.1, 1, 5, 10, 50, 100, 500)),
                            radial = list(cost = c(0.1, 1, 5, 10, 50, 100, 500),
                                         gamma = 2^{-2:2}),
                            polynomial = list(cost = c(0.1, 1, 5, 10, 50, 100, 500),
                                             degree = 2:4))
 }
}

if(length(grids) == 0) stop("'classifier' does not need any tuning \n")
innerlength <- unlist(lapply(grids, length))
if(any(innerlength == 0)) stop("Invalid grids specified \n")

hypergrid <- expand.grid(grids)
hypergrid <- data.frame(apply(hypergrid, 2, sort))

tunereslist <- vector(mode="list", length=nrow(learnmatrix))

if(missing(genesel)){
 for(i in 1:nrow(learnmatrix)){
  if(trace) cat("tuning iteration", i, "\n")
  Xi <- X[learnmatrix[i,],,drop=FALSE]
  yi <- y[learnmatrix[i,]]
  lsi <- GenerateLearningsets(y=yi, method="CV", fold=fold, strat=strat)
  perf <- double(nrow(hypergrid))
  for(k in 1:nrow(hypergrid)){
   classifk <- do.call("classification", args=c(list(X=Xi, y=yi, learningsets=lsi, trace = FALSE,
                                                 classifier = classifier), as.list(data.frame(hypergrid[k,,drop=FALSE])), ll))
   evalk <- evaluation(classifk, scheme = "iterationwise")
   perf[k] <- mean(evalk@score)
 }
  tunereslist[[i]] <- perf
 }
}

else{

 ranks <- genesel@rankings
 imps <- genesel@importance
if(is.element(genesel@method, c("lasso", "elasticnet", "boosting"))){
  if(length(ranks) > 1){
  for(i in 1:nrow(learnmatrix)){
  if(trace) cat("tuning iteration", i, "\n")
  seli <- c()
   for(j in 1:length(ranks)){
    rankj <- ranks[[j]][i,]
    impj <- imps[[j]][i,]
    impj <- impj[impj > 0]
    nbgene <- min(length(impj), nbgene)
    seli <- c(seli, rankj[1:nbgene])
   }
   seli <- unique(seli)
   Xi <- X[learnmatrix[i,],seli,drop=FALSE]
   yi <- y[learnmatrix[i,]]
   lsi <- GenerateLearningsets(y=yi, method="CV", fold=fold, strat=strat)
   perf <- double(nrow(hypergrid))
   for(k in 1:nrow(hypergrid)){
   classifk <- do.call("classification", args=c(list(X=Xi, y=yi, learningsets=lsi, trace = FALSE,
                                                 classifier = classifier), as.list(data.frame(hypergrid[k,,drop=FALSE])), ll))
   evalk <- evaluation(classifk, scheme = "iterationwise")
   perf[k] <- mean(evalk@score)
   }
   tunereslist[[i]] <- perf
  }
 }
 else{
  ranks <- ranks[[1]]
  imps <- imps[[1]]
  for(i in 1:nrow(learnmatrix)){
  if(trace) cat("tuning iteration", i, "\n")
  impi <- imps[i,]
  impi <- impi[impi > 0]
  nbgene <- min(length(impi), nbgene)
  seli <- ranks[i,1:nbgene]
  Xi <- X[learnmatrix[i,],seli,drop=FALSE]
  yi <- y[learnmatrix[i,]]
  lsi <- GenerateLearningsets(y=yi, method="CV", fold=fold, strat=strat)
  perf <- double(nrow(hypergrid))
  for(k in 1:nrow(hypergrid)){
   classifk <- do.call("classification", args=c(list(X=Xi, y=yi, learningsets=lsi, trace = FALSE,
                                                 classifier = classifier), as.list(data.frame(hypergrid[k,,drop=FALSE])), ll))
   evalk <- evaluation(classifk, scheme = "iterationwise")
   perf[k] <- mean(evalk@score)
  }
  tunereslist[[i]] <- perf
  }
 }
}
else{
  if(length(ranks) > 1){
  for(i in 1:nrow(learnmatrix)){
  if(trace) cat("tuning iteration", i, "\n")
  seli <- c()
   for(j in 1:length(ranks)){
    rankj <- ranks[[j]][i,]
    seli <- c(seli,  rankj[1:nbgene])
   }
   seli <- unique(seli)
   Xi <- X[learnmatrix[i,],seli,drop=FALSE]
   yi <- y[learnmatrix[i,]]
   lsi <- GenerateLearningsets(y=yi, method="CV", fold=fold, strat=strat)
   perf <- double(nrow(hypergrid))
   for(k in 1:nrow(hypergrid)){
   classifk <- do.call("classification", args=c(list(X=Xi, y=yi, learningsets=lsi, trace = FALSE,
                                                 classifier = classifier), as.list(data.frame(hypergrid[k,,drop=FALSE])), ll))
   evalk <- evaluation(classifk, scheme = "iterationwise")
   perf[k] <- mean(evalk@score)
   }
   tunereslist[[i]] <- perf
  }
 }
 else{
  ranks <- ranks[[1]]
  for(i in 1:nrow(learnmatrix)){
  if(trace) cat("tuning iteration", i, "\n")
  seli <- ranks[i,1:nbgene]
  Xi <- X[learnmatrix[i,],seli,drop=FALSE]
  yi <- y[learnmatrix[i,]]
  lsi <- GenerateLearningsets(y=yi, method="CV", fold=fold, strat=strat)
  perf <- double(nrow(hypergrid))
  for(k in 1:nrow(hypergrid)){
   classifk <- do.call("classification", args=c(list(X=Xi, y=yi, learningsets=lsi, trace = FALSE,
                                                 classifier = classifier), as.list(data.frame(hypergrid[k,,drop=FALSE])), ll))
   evalk <- evaluation(classifk, scheme = "iterationwise")
   perf[k] <- mean(evalk@score)
   }
   tunereslist[[i]] <- perf
   }
  }
 }
}

 return(new("tuningresult", hypergrid = hypergrid, tuneres = tunereslist,
         method=evalk@method , fold=fold))

})

### X=matrix, y=factor, f=missing

setMethod("tune", signature(X="matrix", y="factor", f="missing"),
          function(X, y, learningsets, genesel, genesellist = list(), nbgene,
                  classifier, fold = 3, strat = FALSE, grids = list(), trace=TRUE, ...){
tune(X, y=as.numeric(y)-1, learningsets=learningsets,
               genesel = genesel, genesellist = genesellist, nbgene = nbgene,
               classifier = classifier, fold = fold, strat = strat,
               grids = grids, trace = trace, ...)
})

### signature X=data.frame, f=formula

setMethod("tune", signature(X="data.frame", y="missing", f="formula"),
          function(X, y, f, learningsets, genesel, genesellist = list(), nbgene,
                  classifier, fold = 3, strat = FALSE, grids = list(), trace=TRUE, ...){
yvar <- all.vars(f)[1]
xvar <- strsplit(as.character(f), split = "~")[[3]]
where <- which(colnames(X) == yvar)
if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
else y <- get(yvar)
if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
f <- as.formula(paste("~", xvar))
X <- model.matrix(f, data=X)[,-1,drop=FALSE]
tune(as.matrix(X), y=y, learningsets=learningsets,
               genesel = genesel, genesellist = genesellist, nbgene = nbgene,
               classifier=classifier, fold = fold, strat = strat,
               grids = grids, trace = trace, ...)})

### X=ExpressionSet, y="character", f="missing"

setMethod("tune", signature(X="ExpressionSet", y="character", f="missing"),
          function(X, y, learningsets, genesel, genesellist = list(), nbgene,
                  classifier, fold = 3, strat = FALSE, grids = list(), trace=TRUE, ...){
          y <- pData(X)[,y]
          X <-  exprs(X)
          if(nrow(X) != length(y)) X <- t(X)
          tune(as.matrix(X), y=y, learningsets=learningsets,
               genesel = genesel, genesellist, nbgene = nbgene, classifier = classifier,
               fold = fold, strat = strat, grids = grids, trace = trace, ...)})