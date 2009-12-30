###**************************************************************************###
### filename: classes.r
### Title: Raw versions of classes.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 18.9.2007

### Brief description:
#
#   Stores the class definitions and convenience methods for them,
#   is frequently updated.
#
### Further comments and notes:
#
#
#
###**************************************************************************###

#+++++++++++ Class: cloutput +++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="cloutput",
        representation(learnind = "numeric", y="numeric", yhat="numeric", 
                        prob="matrix", method="character", mode="character", model="list"))
                       
#+++++++++++ Class: varseloutput ++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="varseloutput",  representation(varsel="numeric"))

### remark: return type of any "filter" method (s. filter. r)
                       
#+++++++++++ sub-Class: clvarseloutput ++++++++++++++++++++++++++++++++++++++#

setClass(Class="clvarseloutput", contains=c("cloutput", "varseloutput"))

### e.g. ElasticNet, Lasso, compBoosting, RandomForest

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

setMethod("show", signature="cloutput", function(object){
          methodprint <- switch(object@method, 
                                compBoost = "Componentwise Boosting",
                                DLDA = "diagonal discriminant analysis",
                                ElasticNet = "Elastic Net",
                                FDA = "Fisher's linear discriminant",
                                flexDA = "flexible discriminant analysis" ,
                                gbm = "Tree-based boosting",
                                knn = "k nearest neighbours",
                                Lasso = "Lasso",
                                LDA = "linear discriminant analysis",
                                nnet = "Feed-forward neural network",
                                pknn = "probabilistic k nearest neighbours",
                                plr = "penalized logistic regression",
                                pls_lda = "partial least squares + lda",
                                pls_lr = "partial least squares + logistic regression",
                                pls_rf = "partial least squares + random forests",
                                pnn = "probabilistic neural networks",
                                QDA = "quadratic discriminant analysis",
                                rf = "random forest",
                                scDA = "shrunken centroids discriminant analysis",
                                shrinkLDA = "shrinkage linear discriminant analysis",
                                svm = "support vector machine")
          if(is.null(methodprint)) methodprint <- object@method                        
          cat(object@mode, " Classification with ", methodprint, "\n", sep="")
          cat("number of predictions: ", length(object@y), "\n", sep="")
          })
          
setGeneric("ftable")
          
setMethod("ftable", signature="cloutput", function(x){
          dims <- sort(unique(c(x@y, x@yhat))) 
          ftab <- matrix(0, nrow=length(dims), ncol=length(dims))
          for(i in seq(along=dims)){
           for(j in seq(along=dims)){
           ftab[i,j] <- length(intersect(which(x@y==dims[i]), which(x@yhat == dims[j])))
           }
          } 
          ftab <- as.table(ftab)
          attr(ftab, "dimnames") <- list(true = as.character(dims), predicted = as.character(dims))
          cat("number of missclassifications: ",sum(ftab)-sum(diag(ftab)), "\n")
          cat("missclassification rate: ", round((sum(ftab)-sum(diag(ftab)))/sum(ftab),3), "\n")
          if(x@mode == "binary"){
          cat("sensitivity:", round(ftab[2,2]/sum(ftab[2,]),3), "\n")
          cat("specificity:", round(ftab[1,1]/sum(ftab[1,]),3), "\n")
          }
          print(ftab)
          cat("\n")})
          
setGeneric("plot")
          
setMethod("plot", signature(x="cloutput", y="missing"), function(x, main = ""){
           if(any(is.na(x@prob)))
           stop("Probability plot cannot be performed because classification method does not provide estimation of class probabilities \n")
           plotprob(x@prob, x@y, main = main)
          })
          
setGeneric("roc", function(object,...) standardGeneric("roc"))

setMethod("roc", signature(object="cloutput"), function(object, plot = TRUE, ...){
           pr <- object@prob
           if(any(is.na(pr))) stop("compuation of the ROC curve requires probabilities \n")
           if(object@mode == "multiclass") stop("computation of the ROC curve only possible for binary classification \n")
           if(nrow(pr) < 2) stop("too few observations to compute an ROC curve \n")
           if(length(unique(object@y)) == 1) stop("All test observations are from the same class \n")
           ROCinternal(test = pr[,2], object@y, plot = plot, ...) })
          
#+++++++++++ Class:  learningsets +++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="learningsets", 
         representation(learnmatrix="matrix", method="character", ntrain="numeric",
                        iter="numeric")
         )
         
setMethod("show", signature="learningsets", function(object){
          cat("learningset mode: ", object@method, "\n")
          cat("number of learningsets: ", object@iter, "\n")
          cat("(maximum) number of observations per learning set: ", object@ntrain, "\n") 
          })
          
#+++++++++++ Class:  genesel  +++++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="genesel", representation(rankings = "list", importance = "list",
                                        method = "character", scheme = "character"))
                                        
setMethod("show", signature="genesel", function(object){
          cat("gene selection performed with '", object@method,"'", "\n", sep="")
          cat("scheme used :'", object@scheme,"'", "\n", sep="")
          if(length(object@rankings) == 1) ngenes <- ncol(object@rankings[[1]])
          else ngenes <- unlist(lapply(object@rankings, ncol))[1]
          cat("number of genes: ", ngenes, "\n")
          if(length(object@rankings) == 1) niter <- nrow(object@rankings[[1]])
          else niter <- unlist(lapply(object@rankings, nrow))[1]
          cat("number of different learningsets: ", niter, "\n")
          })
                                        
setGeneric("toplist", function(object, k=10, iter = 1, show = TRUE, ...) standardGeneric("toplist"))


setMethod("toplist", signature = "genesel", function(object, k=10, iter=1, show=TRUE,...){
          imp <- object@importance
          R <- object@rankings
          topind <- 1:k
          iter <- iter[1]
          if(length(imp) == 1){
          if(k < 0 | k > ncol(R[[1]])) stop("Invalid value for k. \n") 
          ret <- data.frame(index=R[[1]][iter,topind], importance=imp[[1]][iter, topind])
          rownames(ret) <- NULL
          }
          else{
          ll <- length(imp)
          ret <- vector(mode="list")
          nam <- names(imp)
          for(i in 1:ll){
           tempimp <- imp[[i]]
           tempR <- R[[i]]
           if(k < 0 | k > ncol(tempR)) stop("Invalid value for k. \n") 
           ret[[nam[i]]] <- data.frame(index=tempR[iter,topind], importance=tempimp[iter, topind])
           rownames(ret[[nam[i]]]) <- NULL  
           }     
          }
          if(show){
            cat("top ", k, " genes for iteration ", iter, "\n \n")
            print(ret)
            }    
          return(invisible(ret))})
          
          
setMethod("plot", signature(x = "genesel", y = "missing"), function(x, top=10, iter=1, ...){
           imp <- x@importance
           R <- x@rankings
           topind <- 1:top
           iter <- iter[1] 
           dotsCall <- substitute(list(...))
           ll <- eval(dotsCall)
           if(!hasArg(xlab)) ll$xlab <- "gene index" 
           if(!hasArg(ylab)) ll$ylab <- "relative variable Importance"
           if(!hasArg(main)) ll$main <- paste("Variable importance plot for iteration", substitute(iter))
           if(!hasArg(cex.lab)) ll$cex.lab <- 1.5
           if(length(imp) == 1){
           if(top < 0 | top > ncol(R[[1]])) stop("Invalid value for 'top'. \n") 
           imp <- imp[[1]][iter,]
           tempimp <- (imp/sum(imp))[topind]
           plotind <- which(tempimp > 0)
           if(length(plotind) == 0){
            warning("No variable with importance > 0 found. \n")
            next
           }
           tempimp <- (tempimp/sum(tempimp))[intersect(topind, plotind)]
           ll$height <- tempimp
           ll$names.arg <- rep("", length(tempimp))
           deltay <- max(ll$height)/top
           if(!hasArg(ylim)) ll$ylim <- c(0, max(ll$height)+deltay*2)
           bb <- do.call(barplot, args=ll)
           for(i in 1:top){
           chars <- as.character(R[[1]][iter,topind]) 
           characterplot(chars[i], bb[i], ll$height[i], deltax=3/top, deltay=deltay, cex=15/top)}
          }
          else{
          limp <- length(imp)
          nam <- names(imp)
          #plot.new() 
          ask <- ((prod(par("mfcol"))) == 1 && dev.interactive())
          #dev.off()
          opar <- par(ask=ask)
          on.exit(par(opar))
          title <- paste("Variable importance plot for iteration", substitute(iter))
          for(i in 1:limp){
           tempR <- R[[i]]
           if(top < 0 | top > ncol(tempR)) stop("Invalid value for k. \n")
           tempR <- tempR[iter, topind]
           tempimp <- imp[[i]]
           tempimp <- tempimp[iter,]
           plotind <- which(tempimp > 0)
           if(length(plotind) == 0){
            warning("No variable with importance > 0 found. \n")
            next
           }
           tempimp <- (tempimp/sum(tempimp))[intersect(topind, plotind)]
           ll$height <- tempimp
           ll$names.arg <- rep("", length(tempimp))
           deltay <- max(ll$height)/top
           if(!hasArg(ylim)) ll$ylim <- c(0, max(ll$height)+deltay*2)
           if(!hasArg(main)) ll$main <- paste(title, nam[[i]], sep=", ")
           bb <- do.call(barplot, args=ll)
           for(i in 1:top){
           chars <- as.character(tempR) 
           characterplot(chars[i], bb[i], ll$height[i], deltax=4/top, deltay=deltay, cex=15/top)
           }
           }}})
           
#+++++++++++ Class:  evaloutput  ++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="evaloutput", representation(score = "numeric", measure = "character", 
         scheme = "character", method = "character"))
                                        
setMethod("show", signature="evaloutput", function(object){
          cat("evaluated method: '", object@method,"'", "\n", sep="")
          cat("scheme used :'", object@scheme,"'", "\n", sep="")
          cat("peformance measure: '", object@measure, "'", "\n", sep="")
          score <- object@score
          score <- score[!is.na(score)]
          mu <- round(mean(score), digits=3)
          se <- round(sd(score)/sqrt(length(score)), digits=3)
          if(object@scheme != "classwise"){
          cat("mean peformance is", paste(mu), "\n")
          cat("with a standard error of", paste(se), "\n")
          }
          else{
           cat("class-wise performance: \n")
           print(score)
          }})
          
setGeneric("summary", function(object, ...) standardGeneric("summary"))
          
setMethod("summary", signature="evaloutput", function(object){
          if(object@scheme == "classwise") 
          stop("'summary' not meaningful if 'scheme = classwise' \n")
          cat("evaluated method: '", object@method,"'", "\n", sep="")
          cat("scheme used :'", object@scheme,"'", "\n", sep="")
          cat("peformance measure: '", object@measure, "'", "\n", sep="")
          cat("five-point-summary: \n")
          summary(object@score)
          })

setGeneric("boxplot")
          
setMethod("boxplot", signature=(x="evaloutput"), function(x, ...){
            if(x@scheme == "classwise")
             stop("'boxplot' not meaningful if 'scheme = classwise' \n")
            dotsCall <- substitute(list(...))
            ll <- eval(dotsCall)
            if(!hasArg(main)) 
            ll$main=paste(x@method, ": ", x@measure,", ", x@scheme, sep="")
            ll$x <- x@score
            do.call(boxplot, args=ll)
          })

###            
           
setGeneric("obsinfo", function(object, ...) standardGeneric("obsinfo"))
          
setMethod("obsinfo", signature="evaloutput", function(object, threshold=1, show=TRUE){
          if(object@scheme != "observationwise")
          stop("For observationswise misclassification information,
                'scheme' must be 'observationwise' \n")
          if(!is.element(object@measure, c("misclassification", "average probability", "brier score")))
          stop("For observationwise misclassification information, method
                must be one of 'misclassification', 'average probability', 'brier score' \n")
          cat("evaluated method: '", object@method,"'", "\n", sep="")
          cat("scheme used :'", object@scheme,"'", "\n", sep="")
          cat("peformance measure: '", object@measure, "'", "\n \n \n", sep="")
          score <- object@score
          if(object@measure == "average probability") ind <- which(score <= threshold)
          else ind <- which(score >= threshold)
          cat("observations consistently misclassified: \n \n")
          ret <- data.frame(ind, round(score[ind], digits=3))
          colname2 <- ifelse(object@measure == "average probability", object@measure, 
                              paste("average", object@measure))
          colnames(ret) <- c("index", colname2)
          rownames(ret) <- NULL
          indna <- which(is.na(score))
          if(show){
          print(ret)
          cat("\n \n")
          if(length(indna) > 0){
          cat("observations not classified at all: \n \n")
          cat(indna, sep=", ")
          cat("\n \n")
          }
          }
          return(invisible(list(misclassification=ret, notclassified=indna)))
          })

#+++++++++++ Class: tuningresult ++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="tuningresult", representation(hypergrid = "data.frame",
         tuneres = "list", method = "character", fold = "numeric")) 
           
                                        
setMethod("show", signature="tuningresult", function(object){
          cat("tuning for '", object@method,"'", "\n", sep="")
          cat("hyperparameters tuned: \n")
          cat(colnames(object@hypergrid), sep=",")
          cat("\n")
          cat("CV folds used:", object@fold, "\n")
          })

setGeneric("best", function(object, ...) standardGeneric("best"))          


setMethod("best", signature="tuningresult", function(object){
           hypergrid <- object@hypergrid
           tuneres <- object@tuneres
           best <- vector(mode = "list", length = length(tuneres))
           for(i in seq(along = tuneres)){
           bestindi <- which.min(tuneres[[i]])
           temp <- as.list(hypergrid[bestindi,,drop=FALSE])
           attr(temp, "out.attrs") <- NULL
           best[[i]] <- temp
           }
           return(best)})
          


setMethod("plot", signature(x="tuningresult", y="missing"), function(x, iter=1, which=NULL, ...)
          {
          dotsCall <- substitute(list(...))
          ll <- eval(dotsCall)
          if(!hasArg(main)) ll$main <- "results of hyperparameter tuning"
          hypergrid <- x@hypergrid
          tuneres <- x@tuneres
          iter <- iter[1]
          if(is.null(hypergrid) & ncol(hypergrid)>2 || length(which) > 2)
          stop("If more than two hyperparameters have been tuned, 
                at most two can be visualized \n")
          if(!!is.null(hypergrid) && is.element(which, colnames(hypergrid)))
          stop("Names in 'which' do not agree with the names of tuned
                hyperparameters \n")
          if(is.null(which) & ncol(hypergrid) <= 2)
           which <- colnames(hypergrid)
          if(length(which) == 1){
           if(!hasArg(xlab)) ll$xlab <- which
           if(!hasArg(ylab)) ll$ylab <- "misclassification"
           if(!hasArg(ylim)) ll$ylim <- c(0,1)
           ll$x <- hypergrid[,which]
           ll$y <- tuneres[[iter]]
           do.call(plot, args=ll)
           do.call(lines, args=ll)
           bestind <- which.min(tuneres[[iter]])
           abline(h=tuneres[[iter]][bestind], lty="dashed")
           abline(v=hypergrid[bestind,which], lty="dashed")
          }
          else{
           if(!hasArg(xlab)) ll$xlab <- which[1]
           if(!hasArg(ylab)) ll$ylab <- which[2]
           ll$x <- unique(hypergrid[,which[1]])
           ll$y <- unique(hypergrid[,which[2]])
           ll$z <- matrix(tuneres[[iter]], nrow=length(ll$x), ncol=length(ll$y))
           do.call(contour, args=ll)
           bestind <- which.min(tuneres[[iter]])
           points(hypergrid[bestind,which[1]], hypergrid[bestind,which[2]], cex=2, pch=3)
          }
          })
          
		  
		  
#+++++++++++ Class:  predoutput  ++++++++++++++++++++++++++++++++++++++++++++++#
setClass(Class='predoutput',representation(Xnew='matrix',yhat='factor',model='list'))

setMethod("show",signature(object="predoutput"),function(object){
			cat('Number of predictions: ',length(object@yhat),'\n',sep='')
			object@yhat})

setGeneric("prediction",function(X.tr,y.tr,X.new,f,classifier,genesel,models=F,nbgene,tuneres,...) standardGeneric("prediction"))

setMethod("prediction", signature(X.tr='matrix',y.tr='ANY',X.new='matrix',f="missing"), function(X.tr,y.tr,X.new,classifier,genesel,models=F,f,nbgene,tuneres,...){
			
			if(missing(genesel))
				X<-rbind(X.tr,X.new)
			
			if(!missing(genesel)){
				
				if(length(genesel@rankings)!=1) stop('GeneSelection object contains more than one variable selection.')  
				if(missing(nbgene)) stop('nbgene must be specified')
				if(!missing(nbgene)) 
					X<-rbind(X.tr,X.new)[,genesel@rankings[[1]][1:nbgene]]
			}
			
			
			y<-c(y.tr,rep(0,nrow(X.new)))
			learnind<-1:nrow(X.tr)
			
			if(!missing(tuneres)){
				
				if(length(best(tuneres))!=1) stop('bad tuneres object (Contains more than one tuning iteration).')
				
				tu<-best(tuneres)
				
			}
			
			
			dotsCall <- substitute(list(...))
			ll <- eval(dotsCall)
			ll$models<-models
			
			if(!missing(tuneres))
				ll<-c(ll,tu[[1]])
			
			clf<- do.call(classifier, args=c(list(X=X, y=y, learnind = learnind), ll))
			if(is.factor(y.tr)){
				yhn<-clf@yhat+1
				lev<-levels(y.tr)
				yhc<-sapply(yhn,FUN=function(z){lev[z]})
				clf@yhat<-factor(yhc)
			}
			if(!is.factor(y.tr))
				clf@yhat<-factor(clf@yhat)
			    names(clf@yhat)<-1:length(clf@yhat)
			new('predoutput',Xnew=X.new,yhat=clf@yhat,model=clf@model)
			
		})
          
		
###signature X.tr='data.frame', X.new='data.frame',y.tr='missing','f=formula

setMethod("prediction", signature(X.tr='data.frame',y.tr='missing',X.new='data.frame',f="formula"), function(X.tr,y.tr,X.new,f,classifier,genesel,models=F,nbgene,tuneres,...){
          
			yvar <- all.vars(f)[1]
			xvar <- strsplit(as.character(f), split = "~")[[3]]
			where <- which(colnames(X.tr) == yvar)
			if(length(where) > 0 ){  y.tr <- X.tr[,where[1]] ; X.tr<- X.tr[,-where[1]]}
			else stop('Response variable has to be part of data.frames X.tr.')#y.tr <- get(yvar)
			where2<-which(colnames(X.new)==yvar)
			if(length(where)>0){X.new<-X.new[,where2[1]]}
			#f(nrow(X.tr) != length(y.tr)) stop("Number of rows of 'X' must agree with length of y \n")
			f <- as.formula(paste("~", xvar))
			X.tr <- model.matrix(f, data=X.tr)[,-1,drop=FALSE]
			X.new <- model.matrix(f, data=X.new)[,-1,drop=FALSE]
			predicition(X.tr=X.tr,y.tr=y.tr,X.new=X.new,classifier=classifier,genesel=genesel,models=models,f=f,nbgene=nbgene,tuneres=tuneres,...)
			
		})

###signature X.tr="ExpressionSet", X.new="ExpressionSet", y.tr="char", f='missing'.


setMethod("prediction",signature(X.tr="ExpressionSet",y.tr="character",X.new="ExpressionSet",f='missing'),function(X.tr,y.tr,X.new,f,classifier,genesel,models=F,nbgene,tuneres,...){
			
			y.tr <- pData(X.tr)[,y.tr]
			X.tr <-  exprs(X.tr)
			X.new <-  exprs(X.new)
			
			if(nrow(X.tr) != length(y.tr)) { X.tr <- t(X.tr);X.new<-t(X.new)}
			
		prediction(X.tr=X.tr,y.tr=y.tr,X.new=X.new,classifier=classifier,genesel=genesel,models=models,nbgene=nbgene,tuneres=tuneres,...)
			})
                                        

          



                        


          

          
