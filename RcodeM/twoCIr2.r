twoCIr2 <-
function(r2,n,p,lambda=0.95,maxit=100,eps=10^-7,see=1)
{
#Two-sided CI for the multiple square correlation coefficient ro^2 given the observed r2
dump("twoCIr2","c:\\Projects\\Mode\\twoCIr2.r")
#install.packages("hypergo")
library(hypergeo)
#install.packages("nleqslv")
library(nleqslv)
#install.packages("MBESS")
library(MBESS)

fr2=function(r2,ro2,n,p) #pdf of r^2: f(x;ro2)	
{
	fr2=gamma((n-1)/2)/gamma((n-p-1)/2)/gamma(p/2)*(1-ro2)^((n-1)/2)*(r2)^((p-2)/2)*(1-r2)^((n-p-3)/2)
	f0=Re(hypergeo((n-1)/2,(n-1)/2,p/2,r2*ro2))
	fr2=fr2*f0 
	return(fr2)
}
dflog1=function(x,ro2,n,p) #1st derivative of the log pdf: (log f(x;ro2))'	
{
	H0=Re(hypergeo((n-1)/2,(n-1)/2,p/2,x*ro2))
	H1=Re(hypergeo((n-1)/2+1,(n-1)/2+1,p/2+1,x*ro2))
	return(-(n-1)/2/(1-ro2)+x/2*(n-1)^2/p*H1/H0)		
}
cdf=function(r2,ro2,n,p) integrate(fr2,ro2=ro2,n=n,p=p,lower=0,upper=r2)$value #cdf of r^2: F(x;ro2)	
F1=function(r2,ro2,n,p)#1st derivative of the cdf: F'(x;ro2)	 
	integrate(function(x,ro2,n,p) dflog1(x,ro2,n,p)*fr2(x,ro2,n,p),ro2=ro2,n=n,p=p,lower=0,upper=r2)$value

short.nls=function(x,r2,n,pp,lambda=0.95)
{
	#x[1]=r2L;x[2]=r2U
	eq1=cdf(r2=r2,ro2=x[1],n=n,p=pp)-cdf(r2=r2,ro2=x[2],n=n,p=pp)-lambda
	eq2=F1(r2=r2,ro2=x[1],n=n,p=pp)-F1(r2=r2,ro2=x[2],n=n,p=pp)
	c(eq1,eq2)
}

unb.nls=function(x,r2,n,pp,lambda=0.95)
{
	#x[1]=r2L;x[2]=r2U
	eq1=cdf(r2=r2,ro2=x[1],n=n,p=pp)-cdf(r2=r2,ro2=x[2],n=n,p=pp)-lambda
	eq2=fr2(r2=r2,ro2=x[1],n=n,p=pp)-fr2(r2=r2,ro2=x[2],n=n,p=pp)
	c(eq1,eq2)
}
if(see)
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,3,1),cex.lab=1.5,cex.main=1.5)
	r2s=Fs=seq(from=0.01,to=.99,length=100)
	for(i in 1:100) Fs[i]=cdf(r2=r2,ro2=r2s[i],n=n,p=p)
	plot(r2s,Fs,type="l",lwd=3,xlab="ro2",ylab="cdf, probability",main=paste("Cdf as function of ro2 with r2.obs =",r2))
	segments(-1,lambda,1,lambda,lty=2);segments(r2,-1,r2,1)
	segments(-1,1-lambda,1,1-lambda,lty=2)
}

CIout=matrix(0,ncol=2,nrow=5)
#====================Lower-sided CI: the true ro^2 < r2U such that F(r2;r2U)=1-lambda
r2U=r2#start with the estimate r2 itself
for(it in 1:maxit)
{
	F=cdf(r2=r2,ro2=r2U,n=n,p=p)
	if(see) {print(c(it,F-1+lambda,r2U));points(r2U,F,pch=2)}
	delta=(F-1+lambda)/F1(r2=r2,ro2=r2U,n=n,p=p)
	if(abs(delta)<eps) break
	r2U=r2U-delta	
}
CIout[1,2]=r2U

#====================Upper-sided CI: the true ro^2 > r2L such that F(r2;r2L)=lambda
pbeta0=pbeta(r2,p/2,(n-p-1)/2) #cdf F(r2;0)
if(pbeta0>lambda) #compute the lower-sided CI
{	
	r2L=r2#start with the estimate r2 itself
	for(it in 1:maxit)
	{
		F=cdf(r2=r2,ro2=r2L,n=n,p=p)
		if(see) {print(c(it,F-lambda,r2L));points(r2L,F)}
		delta=(F-lambda)/F1(r2=r2,ro2=r2L,n=n,p=p)
		if(abs(delta)<eps) break
		r2L=r2L-delta	
	}
}	
else r2L=0
CIout[2,1]=r2L;CIout[2,2]=1

#====================Equal-tail double-sided CI: the true r2L < ro^2 < r2U where F(r2;r2L)=(1+lambda)/2 and F(r2;r2U)=(1-lambda)/2
r2U=r2#start with the estimate r2 itself
for(it in 1:maxit)
{
	F=cdf(r2=r2,ro2=r2U,n=n,p=p)
	if(see) {print(c(it,F-(1-lambda)/2,r2U));points(r2U,F,pch=2,cex=1.5)}
	delta=(F-(1-lambda)/2)/F1(r2=r2,ro2=r2U,n=n,p=p)
	if(abs(delta)<eps) break
	r2U=r2U-delta	
}
CIout[3,2]=r2U
if(pbeta0>(1+lambda)/2) #compute the lower-sided CI
{	
	r2L=r2#start with the estimate r2 itself
	for(it in 1:maxit)
	{
		F=cdf(r2=r2,ro2=r2L,n=n,p=p)
		if(see) {print(c(it,F-(1+lambda)/2,r2L));points(r2L,F,cex=1.5)}
		delta=(F-(1+lambda)/2)/F1(r2=r2,ro2=r2L,n=n,p=p)
		if(abs(delta)<eps) break
		r2L=r2L-delta	
	}
	CIout[3,1]=r2L
}	
else CIout[3,1]=0



CIout[4,]=CIout[3,] #default CI for unbiased CI
if(CIout[2,1]>0) CIout[4,]=nleqslv(x=as.vector(as.numeric(CIout[3,])),fn=unb.nls,r2=r2,n=n,pp=p,method ="Newton",global="dbldog")$x
CIout[4,1]=max(CIout[4,1],0)
#MBESS CI
out=ci.R2(R2=r2,df.1=p,df.2=n-p-1,Random.Regressors=TRUE,conf.level=lambda)
CIout[5,]=c(out$Lower.Conf.Limit.R,out$Upper.Conf.Limit.R)
if(see) 
{
	segments(CIout[4,],rep(-1,2),CIout[4,],rep(1,2),lty=2,col=3)
	segments(CIout[5,],rep(-1,2),CIout[5,],rep(1,2),lty=2,col=4)
	legend("topright",c("Unbiased","MBESS"),lty=2,col=3:4,cex=1.5)
}
CIout=as.data.frame(CIout,row.names=c("Lower-sided","Upper-sided","Equal-tail","Unbiased","MBESS"))
names(CIout)=c("Lower","Upper")
CIout
}
