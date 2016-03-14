###  GenomeVariancePartition


```R



## Installing BGLR from github
 library(devtools)
 install_git('https://github.com/gdlc/BGLR-R')
 library(BGLR)

## Parameters
  rSqThreshold=.1 # use to determine the end of an LD block
  nQTL=10
  h2=.5
  set.seed(123456)
  
## Data
 data(mice)
 X=mice.X
 p=ncol(X)
 n=nrow(X)

## 'LD-Blocks'
sqCor<-rep(NA,p-1)
blocks<-rep(1,p)

x2=X[,1]
for(i in 2:p){
	x1<-x2
	x2=X[,i]
	sqCor[i]<-cor(x1,x2,use='complete.obs')^2
	if(sqCor[i]>rSqThreshold){
		blocks[i]=blocks[i-1]
	}else{
		blocks[i]=blocks[i-1]+1
	}
}
nBlocks=max(blocks)  

nBlocks   ## ~1,300 blocks from ~ 10,000 markers
p/nBlocks ## in average 7.7 markers per block


## See lots of markers perfectly correlated
  plot(sqCor,type='p',cex=.05,col=2)

## One here means the same block, 0 means change of block
 plot(1-diff(blocks),type='o',col=2,cex=.1,xlim=c(500,1500))


## Simulation
 tmp=floor(p/nQTL)
 whichQTL=seq(from=floor(tmp/2),by=tmp,length=nQTL)
 b0=rep(0,p)
 K=sum(apply(X=X[,whichQTL],FUN=var,MARGIN=2))
 b0[whichQTL]=rnorm(n=nQTL, sd=sqrt(h2/K))
 signal=X%*%b0
 error=rnorm(n,sd=sd(signal)*sqrt((1-h2)/h2))
 y=signal+error
 
## Looking at the 
 varU0_total=var(signal)
 varU_blocks<-rep(NA,nBlocks)
 for(i in 1:nBlocks){
 	varU_blocks[i]=var(X[,blocks==i,drop=F]%*%as.vector(b0[blocks==i]))
 }
 plot(varU_blocks,col=4,type='o');abline(v=blocks[whichQTL],lty=2,col=2)
 
 ## Note: varU_total > sum(varU_blocks) this is likely to be due to un-accounted LD between blocks
 
## Fitting BRR with sets of markers
  ## For examples look at: https://github.com/gdlc/BGLR-R/blob/master/inst/md/example_sets.md

 fm=BGLR(y=y,ETA=list(list(X=X,model='BRR_sets',sets=blocks,saveEffects=T)),nIter=12000,burnIn=2000,thin=3,saveAt='BRR_sets')
 
 plot(fm$ETA[[1]]$varB,cex=.5,col=2); abline(v=whichQTL,col=4,lty=2) # looking at the variance of effects picks only a few and misses many, why?
                                                                    # the ones that we are missing are likely long LD blocks
                                                                    # the average squared effect, this is what the variance is, is small but the block as a whole contributes a lot.
 
 B=readBinMat('BRR_sets_ETA_1_b.bin') # this bring to ram the samples collected for effects across iterations
 
 VAR=getVariances(B=B,X=X,sets=blocks,verbose=T)
 
 
 

```
