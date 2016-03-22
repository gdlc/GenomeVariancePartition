###  GenomeVariancePartition


```R



## Installing BGLR from github
 library(devtools)
 install_git('https://github.com/gdlc/BGLR-R')
 library(BGLR)


### Examples of different approaches for estimating window variances

  library(devtools)
  install_git("https://github.com/gdlc/BGLR-r")
  library(BGLR)
  data(mice); X=mice.X; rm(mice.X,mice.pheno,mice.A)
  
  dir.create('windowVariances')
  setwd('windowVariances')

  ## Parameters
   nQTL=10
   h2=.5
   rSqThreshold=.1
  ##
  
  ## Simulation
   p=ncol(X); n=nrow(X)
   qtl.pos=seq(from=100,by=floor(p/nQTL),length=nQTL)
   b0=rep(0,p)
   b0[qtl.pos]<-sqrt(h2)/sqrt(sum(apply(FUN=var,X=X[,qtl.pos],MARGIN=2)))
  
   signal=X%*%b0
   error=rnorm(sd=sqrt(1-h2),n=n)
   y=signal+error
  ##
  
  
  ## Determination of LD blocks using simple squared correlation
   COR=rep(NA,p)
   blocks=rep(1,p)
   x2=X[,1]
   for(i in 2:p){
   		x1=x2
   		x2=X[,i]
   		COR[i]=cor(x1,x2)
   		if(COR[i]^2<rSqThreshold){
   			blocks[i]=blocks[i-1]+1
   		}else{
   			blocks[i]=blocks[i-1]
   		}	
   }
  

  X=scale(X)
  
 
   ## Running BGLR BayesB, saving effects and then calculating variances per window
  
   fmBB=BGLR(y=y,ETA=list(list(X=X,model='BayesB',saveEffects=T)),
             nIter=6000,burnIn=1000,saveAt='BB_')
   
   BB=readBinMat('BB_ETA_1_b.bin')
   VAR.BB=getVariances(B=BB,X=X,sets=blocks,verbose=T)
  
   plot(colMeans(BB^2))
   abline(v=qtl.pos,col=4,lty=2)
   plot(colMeans(VAR.BB[,-ncol(VAR.BB)]),col=2) #posterior mean of the window variances
   abline(v=blocks[qtl.pos],col=4,lty=2)
  
  
  ## Running BGLR BRR-sets, this method assigns the same variance to markers in the same set
  
   fmBRR_sets=BGLR(y=y,ETA=list(list(X=X,model='BRR_sets',sets=blocks ,saveEffects=T)),
                   nIter=6000,burnIn=1000,saveAt='BRR_sets_') 
   BRR_sets=readBinMat('BRR_sets_ETA_1_b.bin')
   VAR.BRRsets=getVariances(B=BB,X=X,sets=blocks,verbose=T)
  
   plot(colMeans(BRR_sets^2))
   abline(v=qtl.pos,col=2,lty=2)
   
   plot(colMeans(VAR.BRRsets[,-ncol(VAR.BRRsets)]),col=2)
     abline(v=blocks[qtl.pos],col=4,lty=2)
  
   

```
