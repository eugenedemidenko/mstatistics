siBETA <-
function(job=1,n=10,alpha.true=2,beta.true=3,lambda=0.95,N=50,ss=4,nSim=10000,maxiter=100,eps=10^-4)
{
dump(c("siBETA","siBETA2_out"),"c:\\Projects\\Mode\\siBETA.r")
t0=Sys.time()
set.seed(ss)
mle=function(a0,b0,Y,maxiter,eps) # ML estimation of alpha and beta by Fisher scoring algorithm
{
	n=length(Y)
	S1=sum(log(Y));S2=sum(log(1-Y))
	H=matrix(ncol=2,nrow=2)	
	for(it in 1:maxiter)
	{
		dlda=n*digamma(a0+b0)-n*digamma(a0)+S1
		dldb=n*digamma(a0+b0)-n*digamma(b0)+S2
		l=n*lgamma(a0+b0)-n*lgamma(a0)-n*lgamma(b0)+(a0-1)*S1+(b0-1)*S2
		#print(c(it,l,a0,b0,dlda,dldb))
		H[1,1]=trigamma(a0+b0)-trigamma(a0)
		H[1,2]=H[2,1]=trigamma(a0+b0)
		H[2,2]=trigamma(a0+b0)-trigamma(b0)
		covab=-solve(H)/n
		delta=covab%*%c(dlda,dldb)
		if(max(abs(delta))<eps) break
		a0=a0+delta[1]
		b0=b0+delta[2]		
	}
	return(list(l,c(a0,b0),covab))# max log-lik, 2-dim vector of ML alpha and beta, cov matrix
}
polDL=function(ro,aML,bML,th,Y,lambda,nSim) # a function of ro to find the distance from the MLE to the boundary of the DL confidence region
{
	a=max(aML+ro*cos(th),.1);b=max(bML+ro*sin(th),.1)
	n=length(Y)
	fY=prod(dbeta(Y,shape1=a,shape2=b))#the density at observed Y 		
	RB=rbeta(n*nSim,shape1=a,shape2=b)
	DSIM=matrix(dbeta(RB,shape1=a,shape2=b),ncol=n)
	fDSIM=apply(DSIM,1,prod)
	return(mean(fDSIM>=fY)-lambda)
}

Y=rbeta(n,shape1=alpha.true,shape2=beta.true)# create beta-distributed observations
S1=sum(log(Y));S2=sum(log(1-Y))
mY=mean(Y);vY=var(Y)
aMM=mY*(mY*(1-mY)/vY-1);bMM=(1-mY)*(mY*(1-mY)/vY-1) #Method of moments
mle.out=mle(a0=aMM,b0=bMM,Y=Y,maxiter=maxiter,eps=eps)
l0=n*lgamma(alpha.true+beta.true)-n*lgamma(alpha.true)-n*lgamma(beta.true)+(alpha.true-1)*S1+(beta.true-1)*S2

if(job==0)#contours of the density as a function of parameters
{
	aseq=seq(from=alpha.true/5,to=3.5*alpha.true,length=N)#grid over parameter window
	bseq=seq(from=beta.true/5,to=3*beta.true,length=N)
	fY=matrix(ncol=N,nrow=N)
	for(i in 1:N)
	for(j in 1:N)
		fY[i,j]=prod(dbeta(Y,shape1=aseq[i],shape2=bseq[j]))#the density at Y 
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	contour(aseq,bseq,fY,lwd=2,xlab="alpha",ylab="beta")
	points(alpha.true,beta.true,pch=16,cex=1.5)
	Maxf=max(fY)
	contour(aseq,bseq,fY,levels=Maxf*0.995,lwd=2,add=T)
	segments(0,0,100,100,lty=2)
	legend("topleft",paste("n=",n,", alpha.true=",alpha.true,", beta.true=",beta.true,sep=""),cex=1.5,bg="gray95")
}
if(job==1) #example of DL testing & p-value. Comparison of ML with DL
{
	pvML=pchisq(2*(mle.out[[1]]-l0),df=2,lower.tail=F)	#p-value using ML
	
	fY=prod(dbeta(Y,shape1=alpha.true,shape2=beta.true))#the density at observed Y 
		
	RB=rbeta(n*nSim,shape1=alpha.true,shape2=beta.true)
	DSIM=matrix(dbeta(RB,shape1=alpha.true,shape2=beta.true),ncol=n)
	fDSIM=apply(DSIM,1,prod)
	
	r01=runif(n*nSim)# generate uniform r.v. on [0,1]^n
	D01=matrix(dbeta(r01,shape1=alpha.true,shape2=beta.true),ncol=n) #evaluate at pdf
	f01=apply(D01,1,prod) #joint pdf values
	kappa0=0
	for(it in 1:maxiter) #Newton's algotithm
	{
		area=mean(f01>=kappa0)
		LHS=mean(fDSIM>kappa0)-lambda
		delta=LHS/area
		if(abs(delta)<eps) break
		kappa0=kappa0+delta
		#print(c(it,kappa0,LHS))
	}
	pvDL=1-mean(fDSIM>=fY) #p-value using DL
	cat("S1=",S1,", S2=",S2,"\n",sep="")
	cat("pdf observed=",fY,", kappa=",kappa0,sep="","\n")
	cat("p-value ML=",pvML,", p-value DL=",pvDL,sep="","\n")	
}
if(job==2) #simulations
{
	covMLE=covDL=covWL=rep(0,nSim)
	for(isim in 1:nSim)
	{
		Y=rbeta(n,shape1=alpha.true,shape2=beta.true)# create beta-distributed observations
		S1=sum(log(Y));S2=sum(log(1-Y))
		mY=mean(Y);vY=var(Y)
		aMM=mY*(mY*(1-mY)/vY-1);bMM=(1-mY)*(mY*(1-mY)/vY-1) #Method of moments
		mle.out=mle(a0=aMM,b0=bMM,Y=Y,maxiter=maxiter,eps=eps)
		l0=n*lgamma(alpha.true+beta.true)-n*lgamma(alpha.true)-n*lgamma(beta.true)+(alpha.true-1)*S1+(beta.true-1)*S2
		if(2*(mle.out[[1]]-l0)<qchisq(lambda,df=2)) covMLE[isim]=1	
		
		dx=c(alpha.true,beta.true)-mle.out[[2]]
		WML=t(dx)%*%solve(mle.out[[3]])%*%dx
		if(WML<lambda) covWL[isim]=1
		
		fY=prod(dbeta(Y,shape1=alpha.true,shape2=beta.true))#the density at observed Y 		
		RB=rbeta(n*nSim,shape1=alpha.true,shape2=beta.true)
		DSIM=matrix(dbeta(RB,shape1=alpha.true,shape2=beta.true),ncol=n)
		fDSIM=apply(DSIM,1,prod)
		if(mean(fDSIM<fY)<lambda) covDL[isim]=1		
		
	}
	print(Sys.time()-t0)
	covpr3=c(mean(covMLE),mean(covWL),mean(covDL))
	return(covpr3)
	#siBETA2_out=siBETA(job=2)
}
if(job==3) #LR, Wald, and DL confidence curves
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	aseq=seq(from=alpha.true*.1,to=2*alpha.true*1.9,length=N)
	bseq=seq(from=beta.true*.1,to=2*beta.true*1.7,length=N)
#ML confidence regions
	LML=WML=matrix(ncol=N,nrow=N)
	for(i in 1:N)
	for(j in 1:N)
	{
		l0=n*lgamma(aseq[i]+bseq[j])-n*lgamma(aseq[i])-n*lgamma(bseq[j])+(aseq[i]-1)*S1+(bseq[j]-1)*S2
		LML[i,j]=2*(mle.out[[1]]-l0) #LR	
		dx=c(aseq[i],bseq[j])-mle.out[[2]]
		WML[i,j]=t(dx)%*%solve(mle.out[[3]])%*%dx #Wald
	}
	contour(aseq,bseq,LML,levels=qchisq(lambda,df=2),lwd=2,xlim=c(0,8),ylim=c(0,8),xlab="alpha",ylab="beta",drawlabels=F)
	contour(aseq,bseq,WML,levels=qchisq(lambda,df=2),add=T,drawlabels=F)
	
#DL confidence region

	ths=seq(from=0,to=2*pi,length=N);nths=length(ths)
	aML=(mle.out[[2]])[1];bML=(mle.out[[2]])[2]	
	points(aML,bML,pch=17,cex=1.5)
	points(alpha.true,beta.true,pch=16,cex=1.5)
	low.ro=0;up.ro=3*sqrt(mle.out[[3]])[1,1]
	x=y=rep(NA,nths)
	for(i in 1:nths)
	{
		ro=uniroot(f=polDL,aML=aML,bML=bML,th=ths[i],Y=Y,lambda=lambda,nSim=nSim,interval=c(low.ro,up.ro))$root
		x[i]=aML+ro*cos(ths[i]);y[i]=bML+ro*sin(ths[i])
		points(x[i],y[i])
		low.ro=ro*.5;up.ro=1.5*ro	
	}	
	lines(x,y)
	L1=paste("LR (",siBETA2_out[1],")",sep="")
	L2=paste("Wald (",siBETA2_out[2],")",sep="")
	L3=paste("DL (",siBETA2_out[3],")",sep="")
	legend("bottomright",c(L1,L2,L3),lwd=c(2,1,1),pch=c(-1,-1,1),cex=1.5,bg="grey95")
	legend("topleft",c("True parameter","MLE"),pch=c(16,17),,cex=1.5,bg="grey95")
}

print(Sys.time()-t0)	
}
siBETA2_out <-
c(0.93175, 0.3849, 0.95025)
