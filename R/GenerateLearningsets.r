### filename: GenerateLearningsets.r
### Title: Function to prepare different learningsets, e.g.
###        - LOOCV
###        - MCCV
###        - Bootstrap
###
###
### Author: M. Slawski, adapted from A-L Boulesteix
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 19.9.2007
#
### Brief description:
#
#  Returns an object of class learningsets.
#
### Further comments and notes:
#
#   BoostrapCV is, though highly interesting, not included because it
#   does not fit into the general design.
#
###**************************************************************************###

GenerateLearningsets <- function(n, y,
method=c("LOOCV", "CV", "MCCV", "bootstrap"),
         fold=NULL, niter=NULL, ntrain=NULL, strat=FALSE)
{
if(!missing(n)){
 if(length(n) != 1 || n < 0) stop("'n' must be a positive integer ! \n")
 n <- as.integer(n)
 if(!is.null(fold) && n <= fold) stop("'n' is too small \n")
 if(!is.null(ntrain) && n  <= ntrain) stop("'n' is too small \n")
}
if(missing(n) & missing(y)) stop("At least one of 'n' or 'y' mus be given \n")
if(!missing(y)) n <- length(y)

method <- match.arg(method, c("LOOCV","CV","MCCV","bootstrap"))
if(!is.element(method, eval(formals(GenerateLearningsets)$method)))
stop("method must be one of 'LOOCV', 'CV', 'MCCV', 'bootstrap' \n")
if(strat & missing(y))
stop("If 'strat=TRUE', 'y' (class memberships) must be given \n")

if (method=="MCCV")
 {
 if (is.null(niter) | is.null(ntrain))
  stop("With the MCCV method, arguments niter and ntrain should be given.")
   if(strat){
   taby <- table(y)
   prop <- taby/sum(taby)
   classize <- roundvector(prop*ntrain, ntrain)
   if(any(classize < 1))
   stop("Generation of learningsets failed, one or several classes are too small. \n")
   indlist <- sapply(names(taby), function(z) which(y==z), simplify = FALSE)
   learnmatrix <- matrix(nrow=niter, ncol=ntrain)
   lower <- cumsum(c(1, classize[-length(classize)]))
   upper <- cumsum(classize)
   for(i in 1:length(indlist))
   learnmatrix[,lower[i]:upper[i]] <- t(replicate(niter, sample(indlist[[i]], classize[i], replace=FALSE)))
   }

   else learnmatrix <- t(replicate(niter, sample(n,ntrain,replace=FALSE)))
 }

if (method=="CV")
 {
 if (is.null(niter)) niter <- 1
  if (is.null(fold))
  stop("With the CV method, argument 'fold' must be given.")
 if(!strat){
 if (fold==n) method<-"LOOCV"
 else
  {
  size <- n/fold
  #if (size < 5) stop("argument 'fold' is too large; The ratio of no. observations/fold should be > 5. \n")
  learnmatrix <- matrix(0, niter*fold, n-ceiling(size))
  size.int <- floor(size)
  size.vector <- rep(size.int, fold)

  if (size.int != size)
   size.vector[1:((size-size.int)*fold)]<-size.vector[1:((size-size.int)*fold)]+1

  group.index<-c()
  for (j in 1:fold) group.index <- c(group.index, rep(j,size.vector[j]))

  for (i in 1:niter)
   {
   group.index<-group.index[sample(n,n,replace=FALSE)]
   for (j in 1:fold)
    {
    whichj <- which(group.index==j)
    learnmatrix[j+(i-1)*fold,1:length(whichj)]<- whichj
    }
   }
   learnmatrix <- learnmatrix[,1:max(size.vector),drop=FALSE]
   learnmatrix <- t(apply(learnmatrix, 1, function(z) setdiff(0:n, z)))
  }
 }
 else{

   taby <- table(y)
   prop <- taby/sum(taby)
   siz <- n-floor(n/fold)
   classize <- roundvector(prop*siz, siz)
   if(any(taby < fold))
   stop("Generation of learningsets failed, one or several classes are smaller than the number of folds. \n")
   indlist <- sapply(names(taby), function(z) which(y==z), simplify = FALSE)
   templist <- vector(mode="list", length=length(indlist))
   #learnmatrix <- matrix(0, niter*fold, siz)
   #lower <- cumsum(c(1, classize[-length(classize)]))
   #upper <- cumsum(classize)       templist[[i]]
   for(i in 1:length(indlist)){
   outp <- do.call("GenerateLearningsets", args=list(n=taby[i], method="CV", niter=niter, fold=fold))@learnmatrix
   templist[[i]] <- t(apply(outp, 1, function(z) ifelse(z == 0, 0, indlist[[i]][z])))
   #learnmatrix[1:fold,lower[i]:upper[i]] <- t())
   #templist <- lapply(templist, function(z) ifelse(z == 0, 0, indlist[[i]][z]))
   }
   #learnmatrix <- templist[[1]]
   #for(i in 2:length(indlist)) learnmatrix <- cbind(learnmatrix, templist[[i]])
   #checkzeros <- learnmatrix == 0
   #mode(checkzeros) <- "numeric"
   #indmatrix <- rowswaps(checkzeros)
   #for(i in 1:ncol(learnmatrix)) learnmatrix[,i] <- learnmatrix[indmatrix[,i], i]

   topass <- lapply(templist, function(z) z[1:fold,,drop=FALSE])
   swaporder <- rowswaps(topass)
   nrep <- 0
   while(nrep < niter-1){
    swaporder <- rbind(swaporder, swaporder[1:fold,,drop=FALSE]+fold*nrep)
    nrep <- nrep+1
    }



   for(i in 1:length(templist))
   templist[[i]] <- templist[[i]][swaporder[,i],]
   learnmatrix <- templist[[1]]
   for(i in 2:length(indlist)) learnmatrix <- cbind(learnmatrix, templist[[i]])
  }
 }

if (method=="LOOCV") learnmatrix <- matrix(rep(1:n, each=n-1), nrow=n)

if (method=="bootstrap")
 {
 if (is.null(niter))
 stop("If 'method=bootstrap', the argument 'niter' must be given. \n")
 if(!strat) learnmatrix <- t(replicate(niter, sample(n,replace=TRUE)))
 else{
   taby <- table(y)
   if(any(taby) < 1)
   stop("Generation of learningsets failed, one or several classes are too small. \n")
   indlist <- sapply(names(taby), function(z) which(y==z), simplify = FALSE)
   learnmatrix <- matrix(nrow=niter, ncol=n)
   lower <- cumsum(c(1,taby[-length(taby)]))
   upper <- cumsum(taby)
   for(i in 1:length(indlist)){
   learnmatrix[,lower[i]:upper[i]] <- t(replicate(niter, sample(indlist[[i]], taby[i], replace=TRUE)))
  }
 }
 }
 if(strat & is.element(method, c("CV","MCCV","bootstrap")))
 method <- paste("stratified", method)
 new("learningsets", learnmatrix=learnmatrix, method=method,
      ntrain=ncol(learnmatrix), iter=nrow(learnmatrix))
}