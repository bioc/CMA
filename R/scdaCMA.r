### filename: scdaCMA.r
### Title: One of many classifiers.
###
### Author: M. Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 18.9.2007
#
### Brief description:
#
#  Returns an object of class cloutput.
#
### Further comments and notes:
#
#
#
###**************************************************************************###

# Corrected "scdaCMA"-function (all changes are marked with "#########"):

setGeneric("scdaCMA", function(X, y, f, learnind, delta = 0.5,models=FALSE, ...)
			standardGeneric("scdaCMA"))


### signature X=matrix, y=numeric, f=missing:
setMethod("scdaCMA", signature(X="matrix", y="numeric", f="missing"),
		function (X, y, f, learnind, delta = 0.5, models = FALSE, ...) 
		{
			
			nrx <- nrow(X)
			ly <- length(y)
			if (nrx != length(y)) 
				stop("Number of rows of 'X' must agree with length of y \n")
			if (missing(learnind)) 
				learnind <- 1:nrx
			if (length(learnind) > nrx) 
				stop("length of 'learnind' must be smaller than the number of observations. \n")
			y <- as.factor(y)
			K <- nlevels(y)
			levels(y) <- 1:K
			if (K > 2) 
				mode <- "multiclass" else mode <- "binary"
			y <- as.numeric(y) - 1
			Ylearn <- y[learnind]
			Xlearn <- X[learnind, ]
			centroids <- matrix(nrow = ncol(X), ncol = K)
			variances <- matrix(nrow = ncol(X), ncol = K)
			mk <- as.numeric(sqrt(1/table(Ylearn) + 1/length(learnind))) ######### Old faulty version:  mk <- as.numeric(sqrt(1/table(Ylearn) - 1/length(learnind)))
			for (k in 1:K) {
				indk <- (Ylearn == (k - 1))
				centroids[, k] <- muk <- colMeans(Xlearn[indk, , drop = FALSE])
				variances[, k] <- colSums((sweep(Xlearn[indk, , drop = FALSE], 2, muk, "-"))^2) 
			}
			variances <- rowSums(variances)/(length(learnind) - K)
			overallmu <- colMeans(Xlearn)
			D <- sweep(sweep(sweep(centroids, 1, overallmu, "-"), 2, mk, "/"), 1, variances^(-0.5), "*") 
			shrunkD <- sign(D) * pmax(abs(D) - delta, 0)
			shrunkcentroids <- sweep(sweep(sweep(shrunkD, 2, mk, "*"), 1, variances^0.5, "*"), 1, overallmu, "+") 
			priors <- as.numeric(-2 * log(table(Ylearn)/sum(table(Ylearn))))
			Xtest <- X[-learnind, , drop = FALSE]
			if (nrow(Xtest) == 0) {
				Xtest <- Xlearn
				y <- Ylearn
			} else y <- y[-learnind]
			dist <- t(apply(Xtest, 1, function(z)   colSums(sweep(sweep(-shrunkcentroids, 1, z, "+")^2, 1, variances, "/")) + priors))
			dist <- sweep(dist, 1, rowMeans(dist), "-")
			prob <- safeexp(-0.5 * dist)
			prob <- sweep(prob, 1, rowSums(prob), "/")
			yhat <- apply(prob, 1, which.max) - 1
			if (models == TRUE) 
				modd <- list(NULL)
			if (models == FALSE) 
				modd <- list(NULL)
			new("cloutput", yhat = yhat, y = y, learnind = learnind, 
					prob = prob, method = "scDA", mode = mode, model = modd)
		})




### signature X=matrix, y=factor, f=missing:

setMethod("scdaCMA", signature(X="matrix", y="factor", f="missing"),
		function(X, y, learnind, delta = 0.5, models=FALSE, ...){
			scdaCMA(X, y=as.numeric(y)-1, learnind=learnind, delta = delta, models=models,...)
		})

### signature X=data.frame, f=formula

setMethod("scdaCMA", signature(X="data.frame", y="missing", f="formula"),
		function(X, y, f, learnind, delta = 0.5, models=FALSE,...){
			yvar <- all.vars(f)[1]
			xvar <- strsplit(as.character(f), split = "~")[[3]]
			where <- which(colnames(X) == yvar)
			if(length(where) > 0 ){  y <- X[,where[1]] ; X <- X[,-where[1]]}
			else y <- get(yvar)
			if(nrow(X) != length(y)) stop("Number of rows of 'X' must agree with length of y \n")
			f <- as.formula(paste("~", xvar))
			X <- model.matrix(f, data=X)[,-1,drop=FALSE]
			scdaCMA(as.matrix(X), y=y, learnind=learnind, delta = delta, models=models,...)})


### signature: X=ExpressionSet, y=character.

setMethod("scdaCMA", signature(X="ExpressionSet", y="character", f="missing"),
		function(X, y, learnind, delta = 0.5, models=FALSE,...){
			y <- pData(X)[,y]
			X <-  exprs(X)
			if(nrow(X) != length(y)) X <- t(X)
			scdaCMA(X=X, y=y, learnind=learnind, delta = delta, models=models,...)})

