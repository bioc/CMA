###generic
setGeneric("weighted.mcr",function(classifiers,parameters,nbgenes,sel.method,X,y,portion,niter=100,shrinkage=F)standardGeneric("weighted.mcr"))

setMethod("weighted.mcr",signature(classifiers='character',parameters='character',nbgenes="numeric",sel.method='character',X='matrix',y='numeric'),function(classifiers,parameters,nbgenes,sel.method,X,y,portion,niter=100,shrinkage=F){
			
			###check inputs  
			
			if((length(classifiers)== length(parameters) & length(classifiers)==length(nbgenes))==FALSE)
				stop('The arguments classifiers, parameters and nbgenes do not have the same length.')
			
			if((portion<1 & portion>0)==F)
				stop('Invalid input for parameter portion.')
			
			
			
			
#X, y, sel.method and niter will be checked by 'GenerateLearningsets' or 'classification'
			
			###number of training and test data
			n.tr<-ceiling(length(y)*portion)
			n.ts<-length(y)-n.tr
			
			###generate learningsets
			ls<-GenerateLearningsets(y=y,method='MCCV',ntrain=n.tr,niter=niter)
			
			###variable selection
			gsel<-NULL
			if(sel.method!='none')
				gsel<-GeneSelection(y=y,X=X,learningsets=ls,method=sel.method)
			
			###perform resamplings for all classifiers
#matrix of results
			res.mat<-c()
			for(ind in 1:length(classifiers)){
				if(sel.method!='none')
					cli<-  eval(parse(text=paste('classification(classifier=',classifiers[ind],',y=y,X=X,','nbgene=',nbgenes[ind],',learningsets=ls,genesel=gsel,',parameters[ind],')',sep='')))
				else 
					cli<-  eval(parse(text=paste('classification(classifier=',classifiers[ind],',y=y,X=X,','nbgene=',nbgenes[ind],',learningsets=ls,',parameters[ind],')',sep='')))
				
				evali<-evaluation(cli,measure='misclassification')
				res.mat<-cbind(res.mat,evali@score)
			}
			
			
			###mean mcrs
			mcrs<-colMeans(res.mat)
			
			
			###estimate covarainces and correlations
			require(corpcor)
#cov
			cov<-cov(res.mat)
			vars<-diag(cov)
#cor
			corest<-cor(res.mat)
			
			###handle zero variances and other problems in correlation estimation
#make positive definite
			###cor with constant variable -> =0
#handling correlation matrix estiamtion problems (zero variances) and
#correlations of 1 between two different classifiers
			for(one in 1:nrow(corest))
				for(two in 1:ncol(corest)){
					#correlation with constant is zero
					if(is.nan(corest[one,two]) | is.na(corest[one,two]))
						corest[one,two]<-0
					#correlation of 1 between two different classifiers not sensible
					if((corest[one,two]==1 | is.na(corest[one,two]))& one!=two)
						corest[one,two]<-0.9
					#'selfcorrelation' is 1
					if(is.nan(corest[one,two]) & one==two)
						corest[one,two]<-1
					#correlation of 1 between two different classifiers not sensible
					if(corest[one,two]==1 & one!=two)
						corest[one,two]<-0.9
				}
			
			
			###nadeau variances
			
			rho <- n.ts/(n.ts+n.tr)
			nadcovest<- (1/niter+rho/(1-rho))*vars
			nadcovest[nadcovest==0]<-0.0000001
			nadsigest <- make.positive.definite(diag(sqrt(nadcovest)) %*% corest %*% diag(sqrt(nadcovest)))
			
			
			
			if(shrinkage==F){
				###compute weights
				require(mvtnorm)
				weights<-numeric(length(mcrs))
				for(ind1 in 1:length(mcrs)){
					Trans<-c()
					for(ind2 in (1:length(mcrs))[-ind1]){
						vec<-numeric(length(mcrs))
						vec[ind1]<-1
						vec[ind2]<--1
						Trans<-rbind(Trans,vec)
					}
#T<-solve(Trans)
					T<-Trans
					meanst<-T%*%mcrs
					sigmat<-T%*%nadsigest%*%t(T)
					weights[ind1]<-pmvnorm(lower=rep(-Inf,length(meanst)),upper=rep(0,length(meanst)),mean=as.numeric(meanst),sigma=sigmat)
				}
				
				
				###compute corrected/weighted mcr
				corrected.mcr <- as.numeric(t(weights)%*%mcrs)
			}
			if(shrinkage==TRUE){
				###wmcs
#means
				means<-mcrs
				meansold<-mcrs
#rawmean
				rawm<-mean(means)
				
#simulate mnv
				sim<-rmvnorm(50000,meansold,sigma=nadsigest)
#k*(s_0)
				indmin<-which.min(meansold)
#compute bias
				simmin<-sim[apply(sim,1,which.min)==indmin,]
				bias<-min(meansold)-mean(simmin[,indmin])
				
#compute shrinkage factor
				if(bias<=0)
					bias<-0
				if(rawm-min(meansold)!=0)
					sfbias<-bias/(rawm-min(meansold))
				if(rawm-min(meansold)==0)
					sfbias<-1
				if(sfbias>=1)
					sfbias<-1
				
#shrunken mean
				meansnew<-means-sfbias* (means-rawm)
#for mean estimation and in weighted average
				means<-meansnew
#vector of weights
				weightsnew<-rep(NA,length(means))
				
				
#loop over classifiers
				for(ind1 in 1:length(means)){
#differences-operator-matrix for classifier ind1
					Diffs<-c()
					for(ind2 in (1:length(means))[-ind1]){
						vec<-numeric(length(means))
						vec[ind1]<-1
						vec[ind2]<--1
						Diffs<-rbind(Diffs,vec)
					}
					
#means of differences
					meansd<-Diffs%*%means
#covariance matrix of differences
					sigmad<-Diffs%*%nadsigest%*%t(Diffs)
#compute probability that this classifier yields minimal mcr (under
#normality assumption)
					weightsnew[ind1]<-pmvnorm(lower=rep(-Inf,length(meansd)),upper=rep(0,length(meansd)),mean=as.numeric(meansd),sigma=sigmad)
				}
				weights<-weightsnew
				corrected.mcr<-as.numeric(t(weightsnew)%*%meansnew)
				
			}
			
			###get best method
			if(sel.method=='none')
				nbgenes<-rep('all',length(classifiers))
			bestind<-which.min(mcrs)
			best.method<-paste(classifiers[bestind],' with ',parameters[bestind],' and ',nbgenes[bestind], ' variables',sep='')
			
			
			###uncorrected mcr
			uncorrected.mcr<-mcrs[bestind]
			
			###add names
			names(corrected.mcr)<-names(uncorrected.mcr)<-best.method
			
			###ranges
			range.means<-range(mcrs)
			range.ncv<-c(mean(apply(res.mat,1,min)),mean(apply(res.mat,1,max)))
			ranges<-data.frame(range.means,range.ncv)
			dimnames(ranges)<-list(c('min','max'),c('range.means','range.ncv'))
			
			new('wmcr.result',corrected.mcr=corrected.mcr,best.method=best.method,mcrs=mcrs,weights=weights,cov=nadsigest,uncorrected.mcr=uncorrected.mcr,ranges=ranges,mcr.m=res.mat,shrinkage=shrinkage)
			
		}
)

setMethod("weighted.mcr",signature(classifiers='character',parameters='character',nbgenes="numeric",sel.method='character',X='matrix',y='factor'),function(classifiers,parameters,nbgenes,sel.method,X,y,portion,niter=100,shrinkage=F){
			weighted.mcr(classifiers=classifiers,parameters=parameters,nbgenes=nbgenes,sel.method=sel.method,X=X,y=as.numeric(y)-1,portion=portion,niter=niter,shrinkage=shrinkage)
			
		})

setMethod("weighted.mcr",signature(classifiers='character',parameters='character',nbgenes="missing",sel.method='character',X='matrix',y='factor'),function(classifiers,parameters,nbgenes,sel.method,X,y,portion,niter=100,shrinkage=F){
			if(sel.method!='none')
				stop('Number of selected genes must be specified.')
			
			nbgenes<-rep(0,length(classifiers))
			weighted.mcr(classifiers=classifiers,parameters=parameters,nbgenes=nbgenes,sel.method=sel.method,X=X,y=as.numeric(y)-1,portion=portion,niter=niter,shrinkage=shrinkage)
			
		})
