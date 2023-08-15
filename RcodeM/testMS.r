testMS <-
function(job=1,n=10,mu0=1,sd0=0.5,alpha=0.05,N=100,ss=4,nSim=100000)
{
dump("testMS","c:\\Projects\\Mode\\testMS.r")
#install.packages(DEoptim)
library(DEoptim)
#install.packages(nleqslv)
library(nleqslv)
t0=Sys.time()
set.seed(ss)
lambda=1-alpha
#Local functions

int.sim=function(n,kappa,mu0,sd0,muA,sdA,nSim) #approximates integral via simulations
{
	S1=rnorm(nSim,mean=mu0,sd=sd0/sqrt(n))
	S2=sd0^2*rchisq(nSim,df=n-1)
	fst.sim=sqrt(n)*dnorm(sqrt(n)*(S1-muA)/sdA)*dchisq(S2/sdA^2,df=n-1)/sdA^3
	return(mean(fst.sim >= kappa))
}
int.an=function(mu,sigma,mu.star,sigma.star,kappa.star,n) #computes integral analytically
{
	H=-0.5*log(n)+0.5*log(2*pi)+(n-1)/2*log(2)+lgamma((n-1)/2)
	G=log(kappa.star)+H+n*log(sigma.star)
	n3=n-3
		
	s20=exp(2*G/n3)
	for(it in 1:10)
	{
		As2=n3/2*log(s20)-s20/2/sigma.star^2			
		s2.new=s20+2*(G-As2)/((n-3)/s20-1/sigma.star^2)
		if(abs(s2.new-s20)<0.0001) break
		s20=s2.new
	}
	a=s2.new
			
	s20=2*n3*sigma.star^2*(log(2*n3*sigma.star^2)-2*G/n3-1)
	for(it in 1:10)
	{
		As2=(n-3)/2*log(s20)-s20/2/sigma.star^2			
		s2.new=s20+2*(G-As2)/((n-3)/s20-1/sigma.star^2)
		if(abs(s2.new-s20)<0.000001) break
		s20=s2.new
	}
	b=s2.new	
	fint=function(s2,mu,sigma,mu.star,sigma.star,kappa.star,n,G)
	{
		As2=(n-3)/2*log(s2)-s2/2/sigma.star^2		
		cs2=mu.star-sigma.star*sqrt(2*(As2-G)/n)
		ds2=mu.star+sigma.star*sqrt(2*(As2-G)/n)
		fint.ret=1/sigma^2*(pnorm(sqrt(n)*(ds2-mu)/sigma)-pnorm(sqrt(n)*(cs2-mu)/sigma))*dchisq(s2/sigma^2,df=n-1)
		return(fint.ret)
	}
	res=integrate(fint,mu=mu,sigma=sigma,mu.star=mu.star,sigma.star=sigma.star,kappa.star=kappa.star,n=n,G=G,lower=a,upper=b)$value
	return(res)	
}
int.area=function(kappa0,mu,sigma,n) #computes the area of the integral domain
{
	H=-0.5*log(n)+0.5*log(2*pi)+(n-1)/2*log(2)+lgamma((n-1)/2)
	G=log(kappa0)+H+n*log(sigma)
	
	s20=exp(2*G/(n-3))
	for(it in 1:10)
	{
		As2=(n-3)/2*log(s20)-s20/2/sigma^2			
		s2.new=s20+2*(G-As2)/((n-3)/s20-1/sigma^2)
		if(abs(s2.new-s20)<0.0001) break
		s20=s2.new
	}
	a=s2.new
	n3=n-3
	s20=2*n3*sigma^2*(log(2*n3*sigma^2)-2*G/n3-1)
	for(it in 1:10)
	{
		As2=(n-3)/2*log(s20)-s20/2/sigma^2			
		s2.new=s20+2*(G-As2)/((n-3)/s20-1/sigma^2)
		if(abs(s2.new-s20)<0.000001) break
		s20=s2.new
	}
	b=s2.new
	
	fint=function(s2,sigma,n,G)
	{
		As2=(n-3)/2*log(s2)-s2/2/sigma^2		
		fint.ret=2*sigma*sqrt(2*(As2-G)/n)
		return(fint.ret)
	}
	res=integrate(fint,sigma=sigma,n=n,G=G,lower=a,upper=b)$value
	return(res)	
}
un.root=function(kappa.star,mu,sigma,mu.star,sigma.star,n,lambda) #function to find kappa using uniroot
{
	int=int.an(mu=mu,sigma=sigma,mu.star=mu.star,sigma.star=sigma.star,kappa.star=kappa.star,n=n)
	return(int-lambda)
}
int.an.DE=function(x,mu,sigma,n,lambda) #returns the 3d vector of abs(equations) to get 0 for DEoptim for unbiased test  
{
	kappa.star=x[1];sigma.star=x[2];mu.star=x[3]
	H=-0.5*log(n)+0.5*log(2*pi)+(n-1)/2*log(2)+lgamma((n-1)/2)
	G=log(kappa.star)+H+n*log(sigma.star)
	
	s20=exp(2*G/(n-3))
	for(it in 1:10)
	{
		As2=(n-3)/2*log(s20)-s20/2/sigma.star^2			
		s2.new=s20+2*(G-As2)/((n-3)/s20-1/sigma.star^2)
		if(abs(s2.new-s20)<0.0001) break
		s20=s2.new
	}
	a=s2.new
	n3=n-3
	s20=2*n3*sigma.star^2*(log(2*n3*sigma.star^2)-2*G/n3-1)
	for(it in 1:10)
	{
		As2=(n-3)/2*log(s20)-s20/2/sigma.star^2			
		s2.new=s20+2*(G-As2)/((n-3)/s20-1/sigma.star^2)
		if(abs(s2.new-s20)<0.000001) break
		s20=s2.new
	}
	b=s2.new
	
	fint=function(s2,mu,sigma,mu.star,sigma.star,kappa.star,n,H,G)
	{
		As2=(n-3)/2*log(s2)-s2/2/sigma.star^2		
		cs2=mu.star-sigma.star*sqrt(2*(As2-G)/n)
		ds2=mu.star+sigma.star*sqrt(2*(As2-G)/n)
		fint.ret=1/sigma^2*(pnorm(sqrt(n)*(ds2-mu)/sigma)-pnorm(sqrt(n)*(cs2-mu)/sigma))*dchisq(s2/sigma^2,df=n-1)
		return(fint.ret)
	}
	res=integrate(fint,mu=mu,sigma=sigma,mu.star=mu.star,sigma.star=sigma.star,kappa.star=kappa.star,n=n,H=H,G=G,lower=a,upper=b)$value
	ret=rep(NA,3)
	ret[1]=res-lambda
	
	fint.sd=function(s2,mu,sigma,mu.star,sigma.star,kappa.star,n,G)
	{
		As2=(n-3)/2*log(s2)-s2/2/sigma.star^2		
		cs2=mu.star-sigma.star*sqrt(2*(As2-G)/n)
		ds2=mu.star+sigma.star*sqrt(2*(As2-G)/n)
		D=-n/sigma+s2/sigma^3
		t1=sqrt(n)*(ds2-mu)/sigma;t2=sqrt(n)*(cs2-mu)/sigma
		F1=pnorm(t1);F2=pnorm(t2)
		f1=dnorm(t1);f2=dnorm(t2)
		fint.ret=((D+1/sigma)*(F1-F2)-(t1*f1-t2*f2)/sigma)*dchisq(s2/sigma^2,df=n-1)
		return(fint.ret)
	}
	ret[2]=integrate(fint.sd,mu=mu,sigma=sigma,mu.star=mu.star,sigma.star=sigma.star,kappa.star=kappa.star,n=n,G=G,lower=a,upper=b)$value
	fint.mu=function(s2,mu,sigma,mu.star,sigma.star,kappa.star,n,G)
	{
		As2=(n-3)/2*log(s2)-s2/2/sigma.star^2		
		cs2=mu.star-sigma.star*sqrt(2*(As2-G)/n)
		ds2=mu.star+sigma.star*sqrt(2*(As2-G)/n)
		t1=sqrt(n)*(ds2-mu)/sigma;t2=sqrt(n)*(cs2-mu)/sigma
		f1=dnorm(t1);f2=dnorm(t2)
		fint.ret=(f2-f1)*dchisq(s2/sigma^2,df=n-1)
		return(fint.ret)
	}
	ret[3]=integrate(fint.mu,mu=mu,sigma=sigma,mu.star=mu.star,sigma.star=sigma.star,kappa.star=kappa.star,n=n,G=G,lower=a,upper=b)$value
	return(sum(abs(ret)))		
}

#DL test
	S1=rnorm(nSim,mean=mu0,sd=sd0/sqrt(n))
	S2=sd0^2*rchisq(nSim,df=n-1)
	kapMAX=max(sqrt(n)*dnorm(sqrt(n)*(S1-mu0)/sd0)*dchisq(S2/sd0^2,df=n-1)/sd0^3)
	kapAPP=(1-lambda)/kapMAX
	kap0s=int.kap=seq(from=0,to=kapMAX,length=N)
	for(ik in 1:N)
		int.kap[ik]=int.sim(n=n,kappa=kap0s[ik],mu0=mu0,sd0=sd0,muA=mu0,sdA=sd0,nSim)
	a=abs(int.kap-lambda)
	kap0=kap00=kap0s[a==min(a)]
		
	#Newton's algorithm via int.area
	for(k in 1:100)
	{
		Fi=int.an(mu=mu0,sigma=sd0,mu.star=mu0,sigma.star=sd0,kappa.star=kap0,n=n) #computes integral
		area=int.area(kappa0=kap0,mu=mu0,sigma=sd0,n=n) #computes integral
		new.kap=kap0+(Fi-lambda)/area
		if(abs(kap0-new.kap)<10^-6) break		
		kap0=new.kap
	}
	#quasi-Newton's algorithm
	kap0M=(1-lambda)*kapMAX
	for(k in 1:100)
	{
		Fi=int.an(mu=mu0,sigma=sd0,mu.star=mu0,sigma.star=sd0,kappa.star=kap0M,n=n) #computes integral
		new.kap=kap0M+(Fi-lambda)*kapMAX
		if(abs(kap0M-new.kap)<10^-6)	break		
		kap0M=new.kap		
	}
	F0=int.an(mu=mu0,sigma=sd0,mu.star=mu0,sigma.star=sd0,kappa.star=kapAPP,n=n)
	F1=int.an(mu=mu0,sigma=sd0,mu.star=mu0,sigma.star=sd0,kappa.star=kap00,n=n)
	F2=int.an(mu=mu0,sigma=sd0,mu.star=mu0,sigma.star=sd0,kappa.star=kap0,n=n)
	kappa.star=uniroot(un.root,mu=mu0,sigma=sd0,mu.star=mu0,sigma.star=sd0,n=n,lambda=lambda,lower=0.001,upper=0.5)$root
	F3=int.an(mu=mu0,sigma=sd0,mu.star=mu0,sigma.star=sd0,kappa.star=kappa.star,n=n)
	F4=int.an(mu=mu0,sigma=sd0,mu.star=mu0,sigma.star=sd0,kappa.star=kap0M,n=n)
	da=as.data.frame(cbind(c(kapAPP,kap00,kap0,kappa.star,kap0M),c(F0,F1,F2,F3,F4)))
	names(da)=c("kappa","integral")
	row.names(da)=c("approximation","simulations","Newton's algorithm","uniroot","quasi-Newton")
	print(da)
	
	muA=seq(from=.3,to=1.7,length=N)
	sdA=seq(from=.1,to=1.1,length=N)
	pow.LR=matrix(ncol=N,nrow=N)
	minPOW=1000
	for(i in 1:N)
	for(j in 1:N)
	{
		pow.LR[i,j]=1-int.an(mu=muA[j],sigma=sdA[i],mu.star=mu0,sigma.star=sd0,kappa.star=kap0,n=n) #computes integral
		if(minPOW>pow.LR[i,j])
		{
			imin=i;jmin=j
			minPOW=pow.LR[i,j]		
		}
	}
	par(mfrow=c(1,1),mar=c(4,4,1,1))
	contour(sdA,muA,pow.LR,levels=c(0.03,0.05,0.1,0.2,0.4,0.6,0.8,0.9),labcex=1)	
	points(sdA[imin],muA[jmin],pch=16,cex=2)
	points(sd0,mu0)
	segments(-1,mu0,10,mu0,lty=2);segments(sd0,-10,sd0,10,lty=2)		
#Unbiased test
	so=DEoptim(fn=int.an.DE,mu=mu0,sigma=sd0,n=n,lambda=lambda,lower=c(.05,.5,.9),upper=c(.07,.7,1.1),control=DEoptim.control(trace=FALSE))
	kappa.star=so$optim$bestmem[1];sigma.star=so$optim$bestmem[2];mu.star=so$optim$bestmem[3]
	
	muA=seq(from=.3,to=1.7,length=N)
	sdA=seq(from=.2,to=1.6,length=N)
	pow.UNB=matrix(ncol=N,nrow=N)
	minUNB=1000
	for(i in 1:N)
	for(j in 1:N)
	{
		pow.UNB[i,j]=1-int.an(mu=muA[j],sigma=sdA[i],mu.star=mu.star,sigma.star=sigma.star,kappa.star=kappa.star,n=n) #computes integral
		if(minUNB>pow.UNB[i,j])
		{
			imin=i;jmin=j
			minUNB=pow.UNB[i,j]		
		}
	}
	contour(sdA,muA,pow.UNB,levels=c(0.075,0.1,0.2,0.4,0.6,0.8,0.9),labcex=1,lty=2,add=T)	
	points(sdA[imin],muA[jmin],pch=17,cex=2)
	mtext(side=1,"s",font=5,cex=2,line=3)
	mtext(side=2,"m",font=5,cex=2,line=2.75)
	legend("topleft",c(paste("DL test (min power = ",round(minPOW,2),")",sep=""),paste("Unbiased test (min power = ",round(minUNB,2),")",sep="")),lty=1:2,pch=16:17,cex=1.5,bg="gray95")
	
	sd.alt=.82;mu.alt=.8
	S1=rnorm(nSim,mean=mu.alt,sd=sd.alt/sqrt(n))
	S2=sd.alt^2*rchisq(nSim,df=n-1)
	
	jf=sqrt(n)/sigma.star^3*dnorm(sqrt(n)*(S1-mu.star)/sigma.star)*dchisq(S2/sigma.star^2,df=n-1)
	pow=mean(jf<=kappa.star)
	points(sd.alt,mu.alt,pch=18,cex=2)
	segments(sd.alt,-100,sd.alt,100,lty=3)
	segments(100,mu.alt,-100,mu.alt,lty=3)
	text(sd.alt,mu.alt+.03,round(pow,2),adj=0)
	print(Sys.time()-t0)
}
