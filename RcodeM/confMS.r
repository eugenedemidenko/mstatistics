confMS <-
function(n=10,mu0=1,sd0=0.5,alpha=0.05,N=100,ss=3,nSim=100000)
{
dump("confMS","c:\\Projects\\Mode\\confMS.r")
#many local functions are the same as in 'testMS'
#install.packages("nleqslv")
library(nleqslv)
t0=Sys.time()
set.seed(ss)
lambda=1-alpha
int.an=function(mu,sigma,mu.star,sigma.star,kappa.star,n) #computes integral
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
int.an.DE.nleqslv=function(x,mu,sigma,n,lambda) #the same as int.an but suiatble for solving using nleqslv package 
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
		#fint.ret=(D*(F1-F2)-(t1*f1-t2*f2)/sigma)*dchisq(s2/sigma^2,df=n-1)/sigma^2
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
		#fint.ret=1/sigma^2*(f2-f1)*dchisq(s2/sigma^2,df=n-1)
		fint.ret=(f2-f1)*dchisq(s2/sigma^2,df=n-1)
		return(fint.ret)
	}
	ret[3]=integrate(fint.mu,mu=mu,sigma=sigma,mu.star=mu.star,sigma.star=sigma.star,kappa.star=kappa.star,n=n,G=G,lower=a,upper=b)$value
	return(ret)		
}

int.area=function(kappa0,mu,sigma,n) #computes integral
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

lambda=1-alpha
	
fS=aalm=INT=LR.as=kappa=kappaUNB=fSUNB=difk=ALTdif=thet.star.sd=thet.star.mu=fSS=matrix(ncol=N,nrow=N)
	
SS1=rnorm(nSim,mean=mu0,sd=sd0/sqrt(n))
SS2=sd0^2*rchisq(nSim,df=n-1)
ci=n*log(sd0^2/SS2*n)+SS2/sd0^2+n*(SS1-mu0)^2/sd0^2-n
cpLRT=mean(ci<qchisq(lambda,df=2))

S1=rnorm(1,mean=mu0,sd=sd0/sqrt(n)) #S1=0.891217 and S2=2.202743 when using ss=3
S2=sd0^2*rchisq(1,df=n-1)
#return(c(S1,S2))

mus=seq(from=.2,to=1.5,length=N)
sds=seq(from=.2,to=1.4,length=N)
	
par(mfrow=c(1,1),mar=c(4,4,1,1))
	
kap0=0.01	
for(i in 1:N)
for(j in seq(from=N,to=1,by=-1))
{
	kappa0=sqrt(n)/sds[i]^3*dnorm(sqrt(n)*(S1-mus[j])/sds[i])*dchisq(S2/sds[i]^2,df=n-1)
	fS[i,j]=1/sds[i]^3*sqrt(n)*dnorm(sqrt(n)*(S1-mus[j])/sds[i])*dchisq(S2/sds[i]^2,df=n-1)
	

	#Newton's algorithm
	for(k in 1:50)
	{
		Fi=try(int.an(mu=mus[j],sigma=sds[i],mu.star=mus[j],sigma.star=sds[i],kappa.star=kap0,n=n))			#computes integral
		if(length(Fi)==0) break			
		area=int.area(kappa0=kap0,mu=mus[j],sigma=sds[i],n=n) #computes integral

		new.kap=kap0+(Fi-lambda)/area
		if(abs(kap0-new.kap)<10^-8) break		
		kap0=new.kap		
	}		
	kappa[i,j]=kap0
	difk[i,j]=Fi-lambda
	ALTdif[i,j]=int.an(mu=mus[j],sigma=sds[i],mu.star=mus[j],sigma.star=sds[i],kappa.star=fS[i,j],n=n)
	INT[i,j]=int.an(mu=mus[j],sigma=sds[i],mu.star=mus[j],sigma.star=sds[i],kappa.star=kappa0,n=n)
	LR.as[i,j]=n*log(sds[i]^2/S2*n)+S2/sds[i]^2+n*(S1-mus[j])^2/sds[i]^2-n
		
	x0=c(kap0,sds[i],mus[j])
	out.nleqslv=nleqslv(x=x0,fn=int.an.DE.nleqslv,mu=mus[j],sigma=sds[i],n=n,lambda=lambda)
	x=out.nleqslv$x
	kappaUNB[i,j]=x[1]		
	fSUNB[i,j]=1/x[2]^3*sqrt(n)*dnorm(sqrt(n)*(S1-x[3])/x[2])*dchisq(S2/x[2]^2,df=n-1)
	thet.star.sd[i,j]=x[2];thet.star.mu[i,j]=x[3]	
}

contour(sds,mus,fS-kappa,levels=0,lty=2,lwd=3,drawlabels=F)#DL boundary
#contour(sds,mus,ALTdif,levels=lambda,col=2,lwd=10,lty=2,drawlabels=F,add=T)#alternative DL boundary
contour(sds,mus,fSUNB-kappaUNB,levels=0,lty=3,lwd=3,col=4,drawlabels=F,add=T)
sdm=contourLines(sds,mus,fSUNB-kappaUNB,levels=0)
sdm.x=sdm[[1]]$x;sdm.y=sdm[[1]]$y
nll=length(sdm.x)
sdm.SO=matrix(0,nrow=nll,ncol=2)
contour(sds,mus,LR.as,levels=qchisq(lambda,df=2),lwd=3,add=T,drawlabels=F)
	
points(sd0,mu0,pch=16,cex=1.5);text(sd0,mu0-.035,"True",font=2)
ss_DL=sqrt(S2/(n-3))
points(ss_DL,S1,pch=18,cex=1.5);text(ss_DL,S1-.035,"DL",font=2)

inn=rep(NA,nSim)

kap=0.01
for(k in 1:50)
{
	Fi=try(int.an(mu=mu0,sigma=sd0,mu.star=mu0,sigma.star=sd0,kappa.star=kap0,n=n))			#computes integral
	if(length(Fi)==0) break			
	area=int.area(kappa0=kap0,mu=mu0,sigma=sd0,n=n) #computes integral

	new.kap=kap0+(Fi-lambda)/area
	if(abs(kap0-new.kap)<10^-5) break		
	kap0=new.kap
}
fS0=1/sd0^3*sqrt(n)*dnorm(sqrt(n)*(SS1-mu0)/sd0)*dchisq(SS2/sd0^2,df=n-1)	
cpBI=mean(fS0>=kap0)
cp.BDL=int.an(mu=mu0,sigma=sd0,mu.star=mu0,sigma.star=sd0,kappa.star=kap0,n=n)
x0=c(qchisq(1-lambda,df=2),sd0,mu0)
out.nleqslv=nleqslv(x=x0,fn=int.an.DE.nleqslv,mu=mu0,sigma=sd0,n=n,lambda=lambda)
x=out.nleqslv$x
	
mM=S1
points(sqrt(S2/n),mM,pch=1,cex=1.5);text(sqrt(S2/n),mM-.035,"MLE/MO",font=2)
	
cp.UDL=int.an(mu=mu0,sigma=sd0,mu.star=x[3],sigma.star=x[2],kappa.star=x[1],n=n)
cA=paste("Asymptotic LRT confidence region(coverage=",round(cpLRT,3),")",sep="")
cB=paste("DL confidence region (coverage=",round(cp.BDL,3),")",sep="")
cU=paste("Unbiased confidence region (coverage=",round(cp.UDL,3),")",sep="")		
text(.7,1.3,"1",cex=1.5)
text(.95,1,"2",cex=1.5)
text(1.18,1.3,"3",cex=1.5)
mtext(side=1,"s",font=5,cex=2,line=3)
mtext(side=2,"m",font=5,cex=2,line=2.75)
legend("bottomright",c(cA,cU,cB),lwd=2,lty=1:3,pch=c("1","2","3"),cex=1.5,bg="gray95")


for(i in 1:N)
fS[i,]=1/sds[i]^3*sqrt(n)*dnorm(sqrt(n)*(S1-mus)/sds[i])*dchisq(S2/sds[i]^2,df=n-1)
#contour(sds,mus,fS-kappa,levels=0,col=3,drawlabels=F,add=T)

M=matrix(ncol=N,nrow=N)
for(i in 1:N)
for(j in 1:N)
{
	fS0=1/sds[i]^3*sqrt(n)*dnorm(sqrt(n)*(S1-mus[j])/sds[i])*dchisq(S2/sds[i]^2,df=n-1)
	M[i,j]=int.an(mu=mu0,sigma=sd0,mu.star=mu0,sigma.star=sd0,kappa.star=fS0,n=n)
}
contour(sds,mus,M,levels=1-alpha,col=2,lwd=3,add=T)
print(Sys.time()-t0)
}
