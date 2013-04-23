###generic
setGeneric("wmc",function(mcr.m,n.tr,n.ts,shrinkage=F)standardGeneric("wmc"))

setMethod("wmc",signature(mcr.m='matrix',n.tr='numeric',n.ts='numeric'),function(mcr.m,n.tr,n.ts,shrinkage=F){
			
#check inputs
			if(any(mcr.m>1))
				stop('Fold errors >1 detected.')
			
			if(any(mcr.m<0))
				stop('Fold errors <0 detected.')
			
			if(n.tr<=1)
				stop('n.tr >>1 required.')
			
			if(n.ts<=1)
				stop('n.ts >1 required.')
			
			if(n.tr%%1!=0 )
				stop('n.tr not an integral value.')
			
			if(n.ts%%1!=0 )
				stop('n.ts not an integral value.')
			
			
			
#load packages
			require(corpcor)
			require(mvtnorm)
			if(shrinkage==F){
#covariance-estimation
				B<-nrow(mcr.m)
				cov<-cov(mcr.m)
				vars <- diag(cov)
#correlation matrix estimate
				corest<-cor(mcr.m)
				
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
				
#nadeau correlation estimation (correlation between folderrors) 
				rho <- n.ts/(n.ts+n.tr)
#nadeau variance estimation
				nadvars<- (1/B+rho/(1-rho))*vars
#handling zero variances
				nadvars[nadvars==0]<-0.0000001
#estimating covariance matrix using new variances estimates and old
#correlation matrix estimate 
				nadsigest <- make.positive.definite(diag(sqrt(nadvars)) %*% corest %*% diag(sqrt(nadvars)))
				
				##compute weights (probability that respective classifier yields minimum mcr)
#vector of means
				means<-colMeans(mcr.m)
#vector of weights
				weights<-numeric(length(means))
				
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
					weights[ind1]<-pmvnorm(lower=rep(-Inf,length(meansd)),upper=rep(0,length(meansd)),mean=as.numeric(meansd),sigma=sigmad)
				}
#compute weighted mean of mean mcrs
				weighted.mcr <- t(weights)%*%means
				
				corrected.mcr<-weighted.mcr
				
			}
			
			if(shrinkage==T){
#covariance-estimation
				B<-nrow(mcr.m)
				cov<-cov(mcr.m)
				vars <- diag(cov)
#correlation matrix estimate
				corest<-cor(mcr.m)
				
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
				
#nadeau correlation estimation (correlation between folderrors) 
				rho <- n.ts/(n.ts+n.tr)
#nadeau variance estimation
				nadvars<- (1/B+rho/(1-rho))*vars
#handling zero variances
				nadvars[nadvars==0]<-0.0000001
#estimating covariance matrix using new variances estimates and old
#correlation matrix estimate 
				nadsigest <- make.positive.definite(diag(sqrt(nadvars)) %*% corest %*% diag(sqrt(nadvars)))
				
				##compute weights (probability that respective classifier yields minimum mcr)
#vector of means
				means<-colMeans(mcr.m)
				
				meansold<-means
				rawm<-mean(means)
				###expected bias nur einmal berechnen
				sim<-rmvnorm(50000,meansold,sigma=nadsigest)
				indmin<-which.min(meansold)
				simmin<-sim[apply(sim,1,which.min)==indmin,]
#ranexpnew<-mean(apply(sim,1,span))
				bias<-min(meansold)-mean(simmin[,indmin])
				
				if(bias<=0)
					bias<-0
				if(rawm-min(meansold)!=0)
					sfbias<-bias/(rawm-min(meansold))
				
				if(rawm-min(meansold)==0)
					sfbias<-1
				
				if(sfbias>=1)
					sfbias<-1
				
				sf<-sfbias
				
				meansnew<-means-sf* (means-rawm)
				
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
				
				wmc.shrink<-t(weightsnew)%*%meansold
				corrected.mcr<-wmc.shrink
				
			}
			
			list(wmc=corrected.mcr,best=which.min(colMeans(mcr.m)),shrinkage=shrinkage)
			
			
		})
