### filename: compare.r
### Title: Convenience function to compare the performance of classifiers.
###
### Author: M. Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 12.10.2007
#
### Brief description:
#
#   - Input is a list of lists as returned by the function 'classification'.
#   - return is a table with rows corresponding to methods
#     and columns corresponding to performance measures
#   - visualization is by boxplots(optional)
#
#
### Further comments and notes:
#
#
#
###**************************************************************************###

setGeneric("compare", function(clresultlist, measure=c("misclassification",
            "sensitivity", "specifity", "average probability", "brier score", "auc"), aggfun = mean, plot = FALSE, ...) standardGeneric("compare"))

setMethod("compare", signature(clresultlist = "list"),
        function(clresultlist, measure = c("misclassification", "sensitivity", "specifity",
                                         "average probability", "brier score", "auc"), aggfun = mean, plot = FALSE, ...){

#if(class(clresultlist) != "list") stop("'clresultlist' must be a list \n")
classes <- unlist(lapply(clresultlist, function(z) unlist(lapply(z, "class"))))
if(any(!extends(classes, "cloutput")))
stop("Incorrect input: 'clresultlist' must be a list whose elements are lists
      of clresultlists of class clouput \n")
lengthes <- unlist(lapply(clresultlist, length))
ll <- unique(lengthes)
if(length(ll) != 1)
stop("All elements of 'clresultlist' must have the same length \n")


col1 <- unlist(lapply(clresultlist, function(z) unique(unlist(lapply(z, slot, "method")))))
if(length(col1) != length(unique(col1)))
stop("No method may occur more than once \n")
perfmatrix <- matrix(nrow = length(col1), ncol=length(measure))
boxplotdata <- vector(mode = "list", length=length(measure))


for(i in seq(along = measure)){
 temp <- matrix(nrow = ll, ncol = length(col1))
 for(j in seq(along = col1)){
 temp[,j] <- evaluation(clresultlist[[j]], measure = measure[i])@score
 perfmatrix[j,i] <- aggfun(temp[,j])
 }
 boxplotdata[[i]] <- temp
}
 colnames(perfmatrix) <- measure
 rownames(perfmatrix) <- col1

 if(plot){
 dotsCall <- substitute(list(...))
 dots <- eval(dotsCall)
 if(!hasArg(names)) dots$names <- col1
 ask <- ((prod(par("mfcol"))) == 1 && dev.interactive())
 opar <- par(ask=ask)
 on.exit(par(opar))
 for(i in seq(along=boxplotdata)){
  if(!hasArg(main)) dots$main <- measure[i]
  dots$x <- data.frame(boxplotdata[[i]])
  do.call("boxplot", args=dots)
 }
}
return(invisible(data.frame(perfmatrix)))
})