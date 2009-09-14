###    Evaluation measures based on the classification results
###
###   filename: evaluation.r
###   Title: Various gene selection methods.
###
###   Author: M. Slawski, adapted from A-L Boulesteix
###   email: <Martin.Slawski@campus.lmu.de>
###   date of creation: 8.10.2007
#
### Brief description:
#
#  Returns an object of class 'evaluationresult'.
#
### Further comments and notes:
#   Also need for tuning purposes.
#
#
###**************************************************************************###

### generic

setGeneric("evaluation", function(clresult, cltrain = NULL, cost =NULL, y=NULL, measure=c("misclassification",
            "sensitivity", "specificity", "average probability", "brier score", "auc", "0.632", "0.632+"),
            scheme = c("iterationwise", "observationwise", "classwise"))
            standardGeneric("evaluation"))

### signature = "list"

setMethod("evaluation", signature(clresult="list"),
          function(clresult, cltrain = NULL, cost = NULL, y = NULL, measure=c("misclassification", "sensitivity", "specificity",
                    "average probability", "brier score", "auc", "0.632", "0.632+"),
                    scheme = c("iterationwise", "observationwise", "classwise")){

if(length(clresult) < 1) stop("'clresult' must contain at least one element \n")
classes <- unlist(lapply(clresult, class))
methods <- unlist(lapply(clresult, slot, name="method"))
if(any(!extends(classes, "cloutput")))
stop("All elements of the list 'clresult' must be of class 'clouput' \n")
if(length(unique(methods)) != 1)
stop("All elements of the list 'clresult' must have been produced by the same method \n")
ll <- length(clresult)
yhatlist <- lapply(clresult, slot, name="yhat")
ylist <- lapply(clresult, slot, name="y")
llyhat <- unlist(lapply(yhatlist, length))
lly <- unlist(lapply(ylist, length))
if(any((llyhat - lly) != 0))
stop("Length of slots 'y' and 'yhat' must agree for each element of 'clresult' \n")

problist <- lapply(clresult, slot, name="prob")

if(any(is.na(problist)) && is.element(measure, c("average probability", "brier score", "auc")))
stop("Classifier to evaluate does not provide estimates for class probabilities. \n
      Choose an appropriate value for 'measure'. \n")

learnindlist <- lapply(clresult, slot, name="learnind")
lulearn <- unlist(lapply(learnindlist, function(z) length(unique(z[z != 0]))))
trainingmode <- identical(lulearn, llyhat)

measure <- match.arg(measure)
if(!is.element(measure, eval(formals(evaluation)$measure)))
stop("Invalid 'measure' specified \n")
scheme <- match.arg(scheme)
if(!is.element(scheme, eval(formals(evaluation)$scheme)))
stop("Invalid 'scheme' specified \n")

if(scheme == "classwise" && is.null(y))
stop("For 'scheme = classwise', argument 'y' must be specified. \n") 


ylvl <- max(unlist(yhatlist))+1
yhatlvl <- max(unlist(ylist))+1

if(is.element(measure, c("sensitivity", "specificity", "auc")) & (ylvl > 2 || yhatlvl > 2))
stop("'sensitivity', 'specificity' or 'auc' are only computed for binary classification \n")

if(measure == "auc" & scheme %in% c("observationwise", "classwise")){
warning("For 'measure=auc' only 'scheme = iterationwise' is possible \n")
 scheme <- "iterationwise"
}

if(measure %in% c("0.632", "0.632+") && scheme == "classwise")
 stop("measures '0.632' or '0.632+'  cannot be combined with 'scheme = classwise' \n") 

if(!is.null(cltrain) && !extends(class(cltrain), "cloutput")) 
stop("'cltrain' must be an object of class 'cloutput' \n")
if(!is.null(cltrain) && cltrain@method != methods[1])
stop("Classfication method for 'cltrain' is incorrect \n")

if(!is.null(cost)){
 dimcost <- try(dim(cost))
 if(inherits(dimcost, "try-error") ||  dimcost[1] != dimcost[2])
 stop("'cost' must be a square matrix \n")
 if(dimcost[1] < max(ylvl, yhatlvl))
 stop("The number of classes exceeds the dimension of the cost matrix \n")
}

if(trainingmode) n <- lulearn
else n <- unique(lulearn + llyhat)

if(length(n) > 1) stop("Error occured while generating the prediction matrix;
                        check that all elements of 'clresult' are based
                        on the same 'learningsets' object  \n")

if(is.element(measure, c("misclassification", "sensitivity", "specificity", "0.632", "0.632+"))){
predmatrix <- matrix(NA, nrow=ll, ncol = n)
 if(measure == "misclassification"){
 if(trainingmode){
  if(is.null(cost)) for (i in 1:ll) predmatrix[i,] <- yhatlist[[i]] != ylist[[i]]
  else for (i in 1:ll) predmatrix[i,] <- apply(cbind(ylist[[i]]+1, yhatlist[[i]]+1), 1, function(z) cost[z[1], z[2]])
  }
  else{
  if(is.null(cost)) for(i in 1:ll) predmatrix[i,-learnindlist[[i]]] <- yhatlist[[i]] != ylist[[i]]
  else for (i in 1:ll) predmatrix[i,-learnindlist[[i]]] <- apply(cbind(ylist[[i]]+1, yhatlist[[i]]+1), 1, function(z) cost[z[1], z[2]]) 
  }
 }
  if(measure == "sensitivity"){
  if(trainingmode)
  for(i in 1:ll) predmatrix[i,] <- ifelse(ylist[[i]] == 0, NA, yhatlist[[i]] == ylist[[i]])
  else
  for(i in 1:ll) predmatrix[i,-learnindlist[[i]]] <- ifelse(ylist[[i]] == 0, NA, yhatlist[[i]] == ylist[[i]])
 }

 if(measure == "specificity"){
  if(trainingmode) for(i in 1:ll) predmatrix[i,] <- ifelse(ylist[[i]] == 1, NA, yhatlist[[i]] == ylist[[i]])
  else
  for(i in 1:ll) predmatrix[i,-learnindlist[[i]]] <- ifelse(ylist[[i]] == 1, NA, yhatlist[[i]] == ylist[[i]])
 }
 
 if(scheme == "iterationwise") score <- rowMeans(predmatrix, na.rm = TRUE)
 if(scheme == "observationwise") score <- colMeans(predmatrix, na.rm = TRUE)
 if(scheme == "classwise"){ 
    splitpredmatrix <- split(predmatrix, as.factor(y))
    score <- unlist(lapply(splitpredmatrix, mean, na.rm =TRUE))
    names(score) <- unique(y)  
 }   

if(measure == "0.632"){
   if(is.null(cltrain)) stop("For the 0.632 estimator, 'cltrain' must be provided \n")
   if(trainingmode) stop("0.632 estimator is not supported if prediction has been performed with the whole learning set \n")
   trainyhat <- cltrain@yhat
   if(length(trainyhat) != n) warning("0.632 estimator has originally been proposed for bootstrap \n")
   for(i in 1:ll)
   predmatrix[i,-learnindlist[[i]]] <- yhatlist[[i]] != ylist[[i]]
   if(scheme == "iterationwise") score <- 0.632*rowMeans(predmatrix, na.rm = TRUE)+(1-0.632)*mean(trainyhat != cltrain@y)
   if(scheme == "observationwise") score <- 0.632*colMeans(predmatrix, na.rm = TRUE)+(1-0.632)*mean(trainyhat != cltrain@y)
 }

 if(measure == "0.632+"){
  if(is.null(cltrain)) stop("For the 0.632+ estimator, 'cltrain' must be provided \n")
  if(trainingmode) stop("0.632+ estimator is not supported if prediction has been performed with the whole learning set \n")
  trainyhat <- cltrain@yhat
   if(length(trainyhat) != n) warning("0.632+ estimator has originally been proposed for bootstrap \n")
   for(i in 1:ll)
   predmatrix[i,-learnindlist[[i]]] <- yhatlist[[i]] != ylist[[i]]
   taby <- table(cltrain@y)
   ordtab <- order(as.numeric(names(taby)))
   taby <- taby[ordtab]
   taby <- taby/sum(taby)
   tabyhat <- c(sum(cltrain@yhat == 0), tabulate(cltrain@yhat, length(taby)-1))
   tabyhat <- tabyhat/sum(tabyhat)
   noinformation <- sum(taby*(1-tabyhat))
   resubs <- mean(trainyhat != cltrain@y)
   if(scheme == "iterationwise") score <- rowMeans(predmatrix, na.rm = TRUE)
   if(scheme == "observationwise") score <- colMeans(predmatrix, na.rm = TRUE)
   score <- pmin(score, noinformation)
   Rhat <- ifelse(score > resubs & noinformation > resubs, (score  - resubs)/(noinformation - resubs), 0)
   w <- 0.632/(1-0.368*Rhat)
   score <- (1-w)*resubs + w*score
 }
}

if(measure == "average probability"){
 lvls <- (0:(ncol(problist[[1]])-1)) 
 predmatrix <- matrix(NA, nrow=ll, ncol = n)
 if(trainingmode){
  for(i in 1:ll){
  G <- model.matrix(~factor(ylist[[i]], levels=lvls)-1)
  predmatrix[i,] <- rowSums(G*problist[[i]])
  }
  }
  else{
   for(i in 1:ll){
   G <- model.matrix(~factor(ylist[[i]], levels=lvls)-1)
   predmatrix[i,-learnindlist[[i]]] <- rowSums(G*problist[[i]])
  }
  }
  if(scheme == "iterationwise") score <- rowMeans(predmatrix, na.rm = TRUE)
  if(scheme == "observationwise") score <- colMeans(predmatrix, na.rm = TRUE)
  if(scheme == "classwise"){ 
    splitpredmatrix <- split(predmatrix, as.factor(y))
    score <- unlist(lapply(splitpredmatrix, mean, na.rm =TRUE))
    names(score) <- unique(y)  
  }   
 }

 if(measure == "brier score"){
 lvls <- (0:(ncol(problist[[1]])-1))
 predmatrix <- matrix(NA, nrow=ll, ncol = n)
 if(trainingmode){
  for(i in 1:ll){
   G <- model.matrix(~factor(ylist[[i]], levels=lvls)-1)
   predmatrix[i,] <- rowSums((G - problist[[i]])^2)
   }
  }
  else{
   for(i in 1:ll){
   G <- model.matrix(~factor(ylist[[i]], levels=lvls)-1)
   predmatrix[i,-learnindlist[[i]]] <- rowSums((G - problist[[i]])^2)
  }
  if(scheme == "iterationwise") score <- rowMeans(predmatrix, na.rm = TRUE)
  if(scheme == "observationwise") score <- colMeans(predmatrix, na.rm = TRUE)
  if(scheme == "classwise"){ 
    splitpredmatrix <- split(predmatrix, as.factor(y))
    score <- unlist(lapply(splitpredmatrix, mean, na.rm =TRUE))
    names(score) <- unique(y)  
  }   
 }
 }

 if(measure == "auc"){
   score <- double(ll)
   for(i in 1:ll){
   probli <- problist[[i]][,2]
   
   
   if(any(is.na(probli)))
   stop(paste('Classifier does not provide class probabilities. AUC cannot be computed. \n',sep=''))
   
if(length(unique(probli))==1){
	warning(paste('Only one distinct probability on learningset ',i,'. AUC set to 0.5 for the corresponding test set. \n',sep='' ))
	score[i]<-0.5
}


   yli <- ylist[[i]]
   if(length(yli) < 2 || length(unique(yli)) < 2)
   stop("measure 'auc' cannot be computed for scheme 'iterationwise';
         too few test observations or all test observations are from
         the same class \n")
if(length(unique(probli))>1)
   score[i] <- ROCinternal(probli, yli, plot=FALSE)
   }
  }
  
return(new("evaloutput", score=score, measure=measure, scheme=scheme,
            method = methods[1]))

})
