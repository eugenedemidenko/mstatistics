binTEST <-
function(job=1,nse=c(20,15),p0se=c(.2,.5),alpha=0.05,N=50,maxit=100,eps=10^-7,nSim=1000000)
{
dump("binTEST","c:\\Projects\\Mode\\binTEST.r")
qLUUNB=function(qLUNB,qUUNB,n,p0,alpha,maxit,eps)
{
	hder=function(n,qa,p0) 
	{
		h=lgamma(n+2)-lgamma(n*qa+1)-lgamma(n*(1-qa)+1)+n*qa*log(p0)+n*(1-qa)*log(1-p0)
		hder=-n*digamma(n*qa+1)+n*digamma(n*(1-qa)+1)+n*log(p0)-n*log(1-p0)
		return(c(h,hder))
	}
	
	H=matrix(ncol=2,nrow=2)
	for(it in 1:maxit)
	{
		rhs1=pbeta(qUUNB,p0*n+1,(1-p0)*n+1)-pbeta(qLUNB,p0*n+1,(1-p0)*n+1)-(1-alpha)
		hL=hder(n=n,qa=qLUNB,p0=p0)
		hU=hder(n=n,qa=qUUNB,p0=p0)
		rhs2=hL[1]-hU[1]
		#print(c(it,qLUNB,qUUNB,rhs1,rhs2))
		H[1,1]=-dbeta(p0,qLUNB*n+1,(1-qLUNB)*n+1)
		H[1,2]=dbeta(p0,qUUNB*n+1,(1-qUUNB)*n+1)
		H[2,1]=hL[2];H[2,2]=-hU[2]
		delta=solve(H)%*%c(rhs1,rhs2)
		if(max(abs(delta))<eps) break
		qLUNB=qLUNB-delta[1];qUUNB=qUUNB-delta[2]
		
	}	
	return(c(qLUNB,qUUNB))
}

if(job==1) #plots power functions of fourt tests
{
	par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),cex.lab=1.5)
	ps=seq(from=0.01,to=0.99,length=N)
	pchs=c(1,2,4,5)
	for(igr in 1:2)
	{
		n=nse[igr];p0=p0se[igr]
		#Discrete/original pbinom cdf
		qLBIN=qbinom(alpha/2,size=n,prob=p0)/n
		qUBIN=qbinom(1-alpha/2,size=n,prob=p0)/n
		Pow=1-pbinom(qUBIN*n,size=n,prob=ps)+pbinom(qLBIN*n,size=n,prob=ps)
		#Clopper-Pearson 	
		qLCP=qbeta(alpha/2,p0*n,(1-p0)*n+1)
		qUCP=qbeta(1-alpha/2,p0*n,(1-p0)*n+1)
		PowCP=1-pbinom(qUCP*n,size=n,prob=ps)+pbinom(qLCP*n,size=n,prob=ps)
		#MET 	
		qLMET=qbeta(alpha/2,p0*n+1,(1-p0)*n+1)
		qUMET=qbeta(1-alpha/2,p0*n+1,(1-p0)*n+1)
		PowMET=1-pbinom(qUMET*n,size=n,prob=ps)+pbinom(qLMET*n,size=n,prob=ps)
		qlu=qLUUNB(qLUNB=qLMET,qUUNB=qUMET,n=n,p0=p0,alpha=alpha,maxit=maxit,eps=eps)
		qLUNB=qlu[1];qUUNB=qlu[2]
		PowUNB=1-pbinom(qUUNB*n,size=n,prob=ps)+pbinom(qLUNB*n,size=n,prob=ps)
	
		matplot(ps,cbind(Pow,PowCP,PowMET,PowUNB),lty=1,type="l",lwd=2,col=1,ylim=c(0,1),ylab="Power, probability",xlab="Alternative probability")
		segments(-1,alpha,2,alpha,lty=2)
		segments(p0,-1,p0,1,lty=2)
		mtext(side=3,paste("n=",n,", p0=",p0,sep=""),cex=1.75,line=.75)
		#Simulations: check the powers using rejection rule	
		pSIM=seq(from=.1,to=.9,by=.1)
		NS=length(pSIM)		
		for(isim in 1:NS)	
		{
			m=rbinom(nSim,size=n,prob=pSIM[isim])
			phat=m/n
			pow.BIN=mean(phat>qUBIN | phat<=qLBIN)
			pow.CP=mean(phat>qUCP | phat<=qLCP)
			pow.MET=mean(phat>qUMET | phat<=qLMET)
			pow.UNB=mean(phat>qUUNB | phat<=qLUNB)
			#pow.UNB=1-mean(phat<=qUUNB & phat>qLUNB & phat>0)
		
			points(rep(pSIM[isim],4),c(pow.BIN,pow.CP,pow.MET,pow.UNB),pch=pchs,cex=1.7)
		}		
		if(igr<2) legend("bottomright",c("BIN","Clopper-Pearson","MET","Unbiased"),lty=1,col=1,pch=pchs,lwd=2,cex=1.4,bg="grey96") 
	}
}
if(job==2) #Example (a)
{
	n=15;m=5;p0=0.5
	p.hat=m/n
	p.bin=2*min(pbinom(m,size=n,prob=p0),pbinom(m,size=n,prob=p0,lower.tail=F))
	p.CP=2*min(pbeta(p.hat,shape1=p0*n,shape2=(1-p0)*n+1),pbeta(p.hat,shape1=p0*n,shape2=(1-p0)*n+1,lower.tail=F))
	p.MET=2*min(pbeta(p.hat,shape1=p0*n+1,shape2=(1-p0)*n+1),pbeta(p.hat,shape1=p0*n+1,shape2=(1-p0)*n+1,lower.tail=F))	
	#Unbiased
	qLUNB=qbeta(alpha/2,shape1=p0*n+1,shape2=(1-p0)*n+1)
	qUUNB=qbeta(1-alpha/2,shape1=p0*n+1,shape2=(1-p0)*n+1)
	hder=function(n,qa,p0) 
	{
		h=lgamma(n+2)-lgamma(n*qa+1)-lgamma(n*(1-qa)+1)+n*qa*log(p0)+n*(1-qa)*log(1-p0)
		hder=-n*digamma(n*qa+1)+n*digamma(n*(1-qa)+1)+n*log(p0)-n*log(1-p0)
		return(c(h,hder))
	}
	H=matrix(ncol=2,nrow=2)
	for(it in 1:maxit)
	{
	
		rhs1=pbeta(qUUNB,p0*n+1,(1-p0)*n+1)-pbeta(qLUNB,p0*n+1,(1-p0)*n+1)-(1-alpha)
		hL=hder(n=n,qa=qLUNB,p0=p0)
		hU=hder(n=n,qa=qUUNB,p0=p0)
		rhs2=hL[1]-hU[1]
		#print(c(it,qLUNB,qUUNB,rhs1,rhs2))
		H[1,1]=-dbeta(p0,qLUNB*n+1,(1-qLUNB)*n+1)
		H[1,2]=dbeta(p0,qUUNB*n+1,(1-qUUNB)*n+1)
		H[2,1]=hL[2];H[2,2]=-hU[2]
		delta=solve(H)%*%c(rhs1,rhs2)
		if(max(abs(delta))<eps) break
		qLUNB=qLUNB-delta[1];qUUNB=qUUNB-delta[2]		
	}
		
	FS=pbinom(m,size=n,prob=p0)
	F0=pbinom(qLUNB*n,size=n,prob=p0)
	if(FS<=F0/alpha) p.UNB=alpha/F0*FS else p.UNB=alpha/(alpha-F0)*FS
	
	pbtR=binom.test(m,n,p0)$p.value
	out=as.data.frame(matrix(cbind(pbtR,p.bin,p.CP,p.MET,p.UNB),ncol=1),row.names=c("binom.test","Binomial","Clopper-Pearson","MET","Unbiased"))
	names(out)="P-value"
	cat("n=",n,", m=",m,", p0=",p0,"\n",sep="")
	print(out)		
}
if(job==3) #Example (b)
{
	n=15;p0=0.5
	p.alt=0.6
	
	out=as.data.frame(matrix(ncol=3,nrow=4),row.names=c("Binomial","Clopper-Pearson","MET","Unbiased"))
	names(out)=c("P-value","Theoretical","Rejection rule")
	
	m=rbinom(nSim,size=n,prob=p.alt)
	phat=m/n
	#Discrete/original binom
	pval.bin=2*pmin(pbinom(m,size=n,prob=p0),pbinom(m,size=n,prob=p0,lower.tail=F))
	pow.p=mean(pval.bin<=alpha)
	qL=qbinom(alpha/2,size=n,prob=p0)/n
	qU=qbinom(1-alpha/2,size=n,prob=p0)/n
	pow.RR=mean(phat>qU | phat<=qL)
	Pow.th=1-pbinom(qU*n,size=n,prob=p.alt)+pbinom(qL*n,size=n,prob=p.alt)
	out[1,1]=pow.p;out[1,2]=Pow.th;out[1,3]=pow.RR
	
	#Clopper-Pearson beta
	p.CP=2*pmin(pbeta(phat,shape1=p0*n,shape2=(1-p0)*n+1),pbeta(phat,shape1=p0*n,shape2=(1-p0)*n+1,lower.tail=F))
	pow.CP=mean(p.CP<=alpha)
	qL=qbeta(alpha/2,shape1=p0*n,shape2=(1-p0)*n+1)
	qU=qbeta(1-alpha/2,shape1=p0*n,shape2=(1-p0)*n+1)
	pow.RR=mean(phat>qU | phat<=qL)
	Pow.th=1-pbinom(qU*n,size=n,prob=p.alt)+pbinom(qL*n,size=n,prob=p.alt)
	out[2,1]=pow.CP;out[2,2]=Pow.th;out[2,3]=pow.RR
	
	#MET beta
	p.MET=2*pmin(pbeta(phat,shape1=p0*n+1,shape2=(1-p0)*n+1),pbeta(phat,shape1=p0*n+1,shape2=(1-p0)*n+1,lower.tail=F))
	pow.MET=mean(p.MET<=alpha)
	qL=qbeta(alpha/2,shape1=p0*n+1,shape2=(1-p0)*n+1)
	qU=qbeta(1-alpha/2,shape1=p0*n+1,shape2=(1-p0)*n+1)
	pow.RR=mean(phat>qU | phat<=qL)
	Pow.th=1-pbinom(qU*n,size=n,prob=p.alt)+pbinom(qL*n,size=n,prob=p.alt)
	out[3,1]=pow.MET;out[3,2]=Pow.th;out[3,3]=pow.RR
	
	
	#Unbiased beta
	qLUNB=qbeta(alpha/2,shape1=p0*n+1,shape2=(1-p0)*n+1)
	qUUNB=qbeta(1-alpha/2,shape1=p0*n+1,shape2=(1-p0)*n+1)
	hder=function(n,qa,p0) 
	{
		h=lgamma(n+2)-lgamma(n*qa+1)-lgamma(n*(1-qa)+1)+n*qa*log(p0)+n*(1-qa)*log(1-p0)
		hder=-n*digamma(n*qa+1)+n*digamma(n*(1-qa)+1)+n*log(p0)-n*log(1-p0)
		return(c(h,hder))
	}
	H=matrix(ncol=2,nrow=2)
	for(it in 1:maxit)
	{
		rhs1=pbeta(qUUNB,p0*n+1,(1-p0)*n+1)-pbeta(qLUNB,p0*n+1,(1-p0)*n+1)-(1-alpha)
		hL=hder(n=n,qa=qLUNB,p0=p0)
		hU=hder(n=n,qa=qUUNB,p0=p0)
		rhs2=hL[1]-hU[1]
		#print(c(it,qLUNB,qUUNB,rhs1,rhs2))
		H[1,1]=-dbeta(p0,qLUNB*n+1,(1-qLUNB)*n+1)
		H[1,2]=dbeta(p0,qUUNB*n+1,(1-qUUNB)*n+1)
		H[2,1]=hL[2];H[2,2]=-hU[2]
		delta=solve(H)%*%c(rhs1,rhs2)
		if(max(abs(delta))<eps) break
		qLUNB=qLUNB-delta[1];qUUNB=qUUNB-delta[2]		
	}
		
	FS=pbeta(phat,shape1=p0*n+1,shape2=(1-p0)*n+1)
	F0=pbeta(qLUNB,shape1=p0*n+1,shape2=(1-p0)*n+1)
	p.UNB=alpha/F0*FS
	p.UNB[FS>F0/alpha]=alpha/(alpha-F0)*(1-FS[FS>F0/alpha])	
	pow.UNB=mean(p.UNB<=alpha)
	pow.RR=mean(phat>qUUNB | phat<=qLUNB)
	Pow.th=1-pbinom(qUUNB*n,size=n,prob=p.alt)+pbinom(qLUNB*n,size=n,prob=p.alt)
	out[4,1]=pow.UNB;out[4,2]=Pow.th;out[4,3]=pow.RR
	cat("n=",n,", p0=",p0,", p.alt=",p.alt,"\n",sep="")	
	print(out)		
}
if(job==4) #Example (c)
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	p0=0.5;p.alt=0.4;pow.alt=0.8;n=150:220
	
	#Bin original
	qL=qbinom(alpha/2,size=n,prob=p0)/n
	qU=qbinom(1-alpha/2,size=n,prob=p0)/n
	Pow.th.BIN=1-pbinom(qU*n,size=n,prob=p.alt)+pbinom(qL*n,size=n,prob=p.alt)
	
	#MET
	qL=qbeta(alpha/2,shape1=p0*n+1,shape2=(1-p0)*n+1)
	qU=qbeta(1-alpha/2,shape1=p0*n+1,shape2=(1-p0)*n+1)
	Pow.th.MET=1-pbinom(qU*n,size=n,prob=p.alt)+pbinom(qL*n,size=n,prob=p.alt)
	
	#Clopper-Pearson
	qL=qbeta(alpha/2,shape1=p0*n,shape2=(1-p0)*n+1)
	qU=qbeta(1-alpha/2,shape1=p0*n,shape2=(1-p0)*n+1)
	Pow.th.CP=1-pbinom(qU*n,size=n,prob=p.alt)+pbinom(qL*n,size=n,prob=p.alt)
	
	#Unbiased
	hder=function(n,qa,p0) 
	{
		h=lgamma(n+2)-lgamma(n*qa+1)-lgamma(n*(1-qa)+1)+n*qa*log(p0)+n*(1-qa)*log(1-p0)
		hder=-n*digamma(n*qa+1)+n*digamma(n*(1-qa)+1)+n*log(p0)-n*log(1-p0)
		return(c(h,hder))
	}
	LN=length(n)
	Pow.th.UNB=rep(NA,LN)
	for(i in 1:LN)
	{
		qLUNB=qbeta(alpha/2,shape1=p0*n[i]+1,shape2=(1-p0)*n[i]+1)
		qUUNB=qbeta(1-alpha/2,shape1=p0*n[i]+1,shape2=(1-p0)*n[i]+1)
		H=matrix(ncol=2,nrow=2)
		for(it in 1:maxit)
		{
			rhs1=pbeta(qUUNB,p0*n[i]+1,(1-p0)*n[i]+1)-pbeta(qLUNB,p0*n[i]+1,(1-p0)*n[i]+1)-(1-alpha)
			hL=hder(n=n[i],qa=qLUNB,p0=p0)
			hU=hder(n=n[i],qa=qUUNB,p0=p0)
			rhs2=hL[1]-hU[1]
			#print(c(it,qLUNB,qUUNB,rhs1,rhs2))
			H[1,1]=-dbeta(p0,qLUNB*n[i]+1,(1-qLUNB)*n[i]+1)
			H[1,2]=dbeta(p0,qUUNB*n[i]+1,(1-qUUNB)*n[i]+1)
			H[2,1]=hL[2];H[2,2]=-hU[2]
			delta=solve(H)%*%c(rhs1,rhs2)
			if(max(abs(delta))<eps) break
			qLUNB=qLUNB-delta[1];qUUNB=qUUNB-delta[2]		
		}
		Pow.th.UNB[i]=1-pbinom(qUUNB*n[i],size=n[i],prob=p.alt)+pbinom(qLUNB*n[i],size=n[i],prob=p.alt)
	}
	X=cbind(Pow.th.BIN,Pow.th.MET,Pow.th.CP,Pow.th.UNB)
	matplot(n,X,type="o",col=1,lty=1,lwd=2,xlab="Sample size, n",ylab="Power")
	nreq=matrix(ncol=2,nrow=4)
	for(i in 1:4)
	{
		xi=X[,i]
		nreq[i,1]=min(n[xi>=pow.alt])
		nreq[i,2]=max(n[xi<=pow.alt])
	}
	text(n[1]+2,0.87,paste("Alternative probability =",p.alt),adj=0,font=3,cex=1.4)
	segments(0,pow.alt,1000,pow.alt,col=2)
	c1=paste("1 Binomial:",nreq[1,1],"-",nreq[1,2],sep="")
	c2=paste("2 MET:",nreq[2,1],"-",nreq[2,2],sep="")
	c3=paste("3 Clopper-Pearson:",nreq[3,1],"-",nreq[3,2],sep="")
	c4=paste("4 Unbiased:",nreq[4,1],"-",nreq[4,2],sep="")
	legend("bottomright",c(c1,c2,c3,c4),cex=1.5,bg="grey96")
}
}
