moder <-
function(job=1,delta=0.1,r=.7,ro=.4,n=10,alpha=.05,maxit=10,eps=10^-5,N=100,Nros=1000,ss=3,nSim=10000) #c.c
{
dump(c("moder","moder42.out"),"c:\\Projects\\Mode\\moder.r")
#install.packages("hypergeo")
set.seed(ss)
t0=Sys.time()
library(hypergeo)


#hypergeo=function(a,b,cc,z)
#{
#	intf=function(x,a,b,cc,z) x^(b-1)*(1-x)^(cc-b-1)/(1-x*z)^a
#	r=integrate(intf,a=a,b=b,cc=cc,z=z,lower=0,upper=1)$value
#	r*gamma(cc)/gamma(b)/gamma(cc-b)
#}
#print(hypergeo(.5,.6,3,-.2))
#return()
#library(pracma)
#hypergeo=function(a,b,cc,z)
#{
#	
#}
lambda=1-alpha
dens.ro=function(r,ro,n) #r density
{
	COF.num=(n-2)*gamma(n-1)*(1-ro^2)^(n/2-.5)*(1-r^2)^(n/2-2)
	COF.den=sqrt(2*pi)*gamma(n-0.5)*(1-ro*r)^(n-1.5)
	fr=COF.num/COF.den*hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5)
	return(Re(fr))
}
lndens.ro=function(ro,r,n) #log r density up to a const
{
	fr=(n-1)/2*log(1-ro^2)-(2*n-3)/2*log(1-ro*r)+log(hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5))
	return(Re(fr))
}
lndens.roDL=function(ro,r,n) #log r density up to a const
{
	fr=(n-4)/2*log(1-r^2)-(2*n-3)/2*log(1-ro*r)+log(hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5))
	return(Re(fr))
}
der1=function(r,ro,n) #derivative of log density
{
	H0=hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5)
	H1=hypergeo(0.5+1,0.5+1,n-0.5+1,ro*r/2+0.5)
	der1=-ro*(n-1)/(1-ro^2)+r*(2*n-3)/2/(1-ro*r)+r/4/(2*n-1)*H1/H0
	return(Re(der1))
}
lndens.der1=function(r,ro,n)
{
	H0=hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5)
	H1=hypergeo(0.5+1,0.5+1,n-0.5+1,ro*r/2+0.5)
	der1=-(n-3)*r/(1-r^2)+(2*n-3)/2*ro/(1-ro*r)+r/4/(2*n-1)*H1/H0
	return(Re(der1))
}
lndens.der1DL=function(r,ro,n)
{
	H0=hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5)
	H1=hypergeo(0.5+1,0.5+1,n-0.5+1,ro*r/2+0.5)
	der1=-(n-4)*r/(1-r^2)+(2*n-3)/2*ro/(1-ro*r)+ro/4/(2*n-1)*H1/H0
	return(Re(der1))
}
der2=function(r,ro,n)
{
	H0=hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5)
	H1=hypergeo(0.5+1,0.5+1,n-0.5+1,ro*r/2+0.5)
	H2=hypergeo(0.5+2,0.5+2,n-0.5+2,ro*r/2+0.5)
	der2=-(1+ro^2)*(n-1)/(1-ro^2)^2+r^2*(2*n-3)/2/(1-ro*r)^2+r^2/16/(2*n-1)*(9/(2*n+1)*H2/H0-1/(2*n-1)*H1^2/H0^2)
	return(Re(der2))
}
der3=function(r,ro,n)
{
	H0=hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5)
	H1=hypergeo(0.5+1,0.5+1,n-0.5+1,ro*r/2+0.5)
	H2=hypergeo(0.5+2,0.5+2,n-0.5+2,ro*r/2+0.5)
	H3=hypergeo(0.5+3,0.5+3,n-0.5+3,ro*r/2+0.5)
	der3=-2*ro*(3+ro^2)*(n-1)/(1-ro^2)^3+r^3*(2*n-3)/(1-ro*r)^3+r^3/64/(2*n-1)*(9/(2*n+1)*(25/(2*n+3)*H3/H0-1/(2*n-1)*H1*H2/H0^2)-1/(2*n-1)*(50/(2*n+3)*H1*H2/H0^2-2/(2*n-1)^2*H1^3/H0^3))
	return(Re(der3))
}
num.ro=function(r,ro,n) (der2(r,ro,n)+(der1(r,ro,n))^2)*dens.ro(r,ro,n)
num.int=function(r,ro,n) integrate(num.ro,ro=ro,n=n,lower=-1,upper=r)$value
den.ro=function(r,ro,n)
{
	d=dens.ro(r,ro,n)
	rr=der3(r,ro,n)*d+3*der1(r,ro,n)*der2(r,ro,n)*d+(der1(r,ro,n))^3*d
	return(rr)
}
dlnd=function(r,ro,n) #derivative of log density * dens
 der1(r,ro,n)*dens.ro(r,ro,n)

den.int=function(r,ro,n) integrate(den.ro,ro=ro,n=n,lower=-1,upper=r)$value


cdf.ro=function(r,ro,n) # cdf
{
	dens.ro=function(x,ro,n)
	{
		COF.num=(n-2)*gamma(n-1)*(1-ro^2)^(n/2-.5)*(1-x^2)^(n/2-2)
		COF.den=sqrt(2*pi)*gamma(n-0.5)*(1-ro*x)^(n-1.5)
		fr=COF.num/COF.den*hypergeo(0.5,0.5,n-0.5,ro*x/2+0.5)
		return(Re(fr))
	}
	cdf.val=try(integrate(f=dens.ro,ro=ro,n=n,lower=-1,upper=r)$value)
	if(!is.numeric(cdf.val)) cdf.val=NA
	return(cdf.val)	
}
cdf.ro.vect=function(r,ro,n)
{
	fint=function(xX)
	{
		dens.ro=function(x,ro,n)
		{
			COF.num=(n-2)*gamma(n-1)*(1-ro^2)^(n/2-.5)*(1-x^2)^(n/2-2)
			COF.den=sqrt(2*pi)*gamma(n-0.5)*(1-ro*x)^(n-1.5)
			fr=COF.num/COF.den*hypergeo(0.5,0.5,n-0.5,ro*x/2+0.5)
			return(Re(fr))
		}
		ret=integrate(dens.ro,n=xX[3],ro=xX[2],lower=-1,upper=xX[1])$value
	}
	X=cbind(r,ro,rep(n,length(r)))
	ret=apply(X,1,FUN=fint)
    return(ret)	
}

densder.ro=function(r,ro,n)
{
#print(length(x))
#print(length(ro))
#print(length(n))
	COF.num=(n-2)*gamma(n-1)*(1-ro^2)^(n/2-.5)*(1-r^2)^(n/2-2)
	COF.den=sqrt(2*pi)*gamma(n-0.5)*(1-ro*r)^(n-1.5)
	f0=hypergeo(0.5,0.5,n-0.5,r*ro/2+0.5)
	fr=Re(COF.num/COF.den*f0)
			
	f1=hypergeo(0.5+1,0.5+1,n-0.5+1,r*ro/2+0.5)
	hr=-ro*(n-1)/(1-ro^2)+r*(2*n-3)/2/(1-ro*r)+r/4/(2*n-1)*Re(f1/f0)
	return(fr*hr)		
}
cdf.dero=function(r,ro,n)  #derivative F'
{
	cdf.val=integrate(f=densder.ro,ro=ro,n=n,lower=-1,upper=r)$value
	return(cdf.val)	
}

covprFisher_fun=function(r,ro,n,alpha=0.05)
{
	Z1a=qnorm(1-alpha/2)
	ex=exp(-2*Z1a/sqrt(n-3))
	rr=(1+r)/(1-r)
	low.ro=(rr*ex-1)/(rr*ex+1)
	ex=exp(2*Z1a/sqrt(n-3))
	up.ro=(rr*ex-1)/(rr*ex+1)
	pr=(ro < up.ro & ro > low.ro)
	return(pr*dens.ro(r=r,ro=ro,n=n))
}
widFisher_fun=function(r,ro,n,alpha=0.05)
{
	Z1a=qnorm(1-alpha/2)
	ex=exp(-2*Z1a/sqrt(n-3))
	rr=(1+r)/(1-r)
	low.ro=(rr*ex-1)/(rr*ex+1)
	ex=exp(2*Z1a/sqrt(n-3))
	up.ro=(rr*ex-1)/(rr*ex+1)
	return((up.ro-low.ro)*dens.ro(r=r,ro=ro,n=n))
}
covprET_fun=function(r,ro,n,alpha=0.05)
{
	Z1a=qnorm(1-alpha/2)
	ex=exp(-2*Z1a/sqrt(n-3))
	rr=(1+r)/(1-r)
	low.ro=(rr*ex-1)/(rr*ex+1)
	ex=exp(2*Z1a/sqrt(n-3))
	up.ro=(rr*ex-1)/(rr*ex+1)
	
	qLET=low.ro;qUET=up.ro
	nr=length(r)
	intr=rep(NA,nr)
	for(i in 1:nr)
	{
		for(it in 1:10)
		{
			delta1=(cdf.ro(r=r[i],qUET[i],n=n)-alpha/2)/integrate(dlnd,ro=qUET[i],n=n,lower=-1,upper=r[i])$value
			delta2=(cdf.ro(r=r[i],qLET[i],n=n)-(1-alpha/2))/integrate(dlnd,ro=qLET[i],n=n,lower=-1,upper=r[i])$value
			qUET[i]=qUET[i]-delta1;qLET[i]=qLET[i]-delta2		
		}
		intr[i]=(ro < qUET[i] & ro > qLET[i])	
	}
	return(intr*dens.ro(r=r,ro=ro,n=n))
}
covprSHORT_fun=function(r,ro,n,alpha=0.05)
{
	Z1a=qnorm(1-alpha/2)
	ex=exp(-2*Z1a/sqrt(n-3))
	rr=(1+r)/(1-r)
	low.ro=(rr*ex-1)/(rr*ex+1)
	ex=exp(2*Z1a/sqrt(n-3))
	up.ro=(rr*ex-1)/(rr*ex+1)
	#Shortest CI
	qL=low.ro;qU=up.ro
	nr=length(r)
	H=matrix(ncol=2,nrow=2)
	intr=rep(NA,nr)
	for(i in 1:nr)
	{
		for(it in 1:30) 
		{
			rh1=cdf.ro(r[i],qL[i],n)-cdf.ro(r[i],qU[i],n)-(1-alpha)
			rh2=cdf.dero(r[i],qL[i],n)-cdf.dero(r[i],qU[i],n)
			H[1,1]=cdf.dero(r[i],qL[i],n);H[1,2]=-cdf.dero(r[i],qU[i],n)
			H[2,1]=num.int(r[i],qL[i],n);H[2,2]=-num.int(r[i],qU[i],n)
			delta=solve(H)%*%c(rh1,rh2)
			q12=c(qL[i],qU[i])-delta			
			if(max(abs(delta))<10^-4) break
			qL[i]=q12[1];qU[i]=q12[2]			
		}
		intr[i]=(ro < qU[i] & ro > qL[i])	
	}	
	return(intr*dens.ro(r=r,ro=ro,n=n))
}
if(job==0) #cdf of r and its approximtions
{
	par(mfrow=c(2,2),mar=c(2,2,2,1),omi=c(.25,.25,0,0))
	ros=c(-.8,-.3,0,.6,.9)
	ns=c(20,10,8,16,30)
	x=cdf_r=cdf_t=cdf_Z=seq(from=-.99,to=0.99,length=N)
	for(i in 1:4)
	{
		ro=ros[i];n=ns[i]		
		X=matrix(rnorm(nSim*n),ncol=n);xm=rowMeans(X)
		X2=rowSums(X^2)
		Z=matrix(rnorm(nSim*n),ncol=n)
		Y=ro*X+Z*sqrt(1-ro^2)
		Y2=rowSums(Y^2);XY=rowSums(X*Y)
		ym=rowMeans(Y)
		r.obs=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))
		r.obs=sort(r.obs)
		for(j in 1:N)
		{
			cdf_r[j]=cdf.ro(r=x[j],ro=ro,n=n)
			xtr=x[j]*sqrt(n-2)/sqrt(1-x[j]^2)
			cdf_t[j]=pt(xtr,df=n-2,ncp=ro*sqrt(n-2)/sqrt(1-ro^2))
			ztr=0.5*log((1+x[j])/(1-x[j]))
			cdf_Z[j]=pnorm(ztr,mean=0.5*log((1+ro)/(1-ro)),sd=1/sqrt(n-3))
		}
		plot(r.obs,(1:nSim)/nSim,xlim=c(-1,1),ylim=c(0,1),type="s",xlab="",ylab="")
		mtext(side=3,paste("n=",n,", ro=",ro,sep=""),cex=1.25,line=0.3)
		lines(x,cdf_r,col=2)
		lines(x,cdf_t,col=3)
		lines(x,cdf_Z,col=4)
		if(i<4) legend("bottomright",c("Simulation-based","Integral over the density","Noncentral t-distribution","Fisher Z-transformation"),col=1:4,lty=1,lwd=2,cex=1,bg="grey96")
		else legend("topleft",c("Simulation-based","Integral over the density","Noncentral t-distribution","Fisher Z-transformation"),col=1:4,lty=1,lwd=2,cex=1,bg="grey96")
	}	
	mtext(side=1,"x",cex=1.5,outer=T)
	mtext(side=2,"cdf",cex=1.5,outer=T)
}
if(job==0.1) #check der1,2,3
{
	ros=lnf=der1.f=intder1=intder2=der2.f=intder3=seq(from=-.9,to=.9,by=.1)
	N=length(ros)
	par(mfrow=c(1,3))
	for(i in 1:N)
	{
		lnf[i]=log(dens.ro(r,ros[i],n))
		intder1[i]=integrate(der1,r=r,n=n,lower=-.999,upper=ros[i])$value	
		der1.f[i]=der1(r=r,ro=ros[i],n=n)
		intder2[i]=integrate(der2,r=r,n=n,lower=-.999,upper=ros[i])$value			
		der2.f[i]=der2(r=r,ro=ros[i],n=n)
		intder3[i]=integrate(der3,r=r,n=n,lower=-.999,upper=ros[i])$value			
	}
	plot(ros,lnf-intder1-mean(lnf-intder1),type="l")
	plot(ros,der1.f-intder2-mean(der1.f-intder2),type="l")
	plot(ros,der2.f-intder3-mean(der2.f-intder3),type="l")
}
if(job==0.2) #check der1,2,3
{
	ros=lnf=der1.f=der2.f=der3.f=seq(from=-.9,to=.9,length=N)
	N=length(ros);N1=N-1
	par(mfrow=c(1,3))
	for(i in 1:N)
	{
		lnf[i]=log(dens.ro(r,ros[i],n))
		der1.f[i]=der1(r=r,ro=ros[i],n=n)
		der2.f[i]=der1(r=r,ro=ros[i],n=n)
		der3.f[i]=der3(r=r,ro=ros[i],n=n)
	}
	h=(max(ros)-min(ros))/N
	matplot(ros[1:N1],cbind(der1.f[1:N1],(lnf[2:N]-lnf[1:N1])/h),type="l")
	matplot(ros[1:N1],cbind(der2.f[1:N1],(der1.f[1:N1]-der1.f[2:N])/h),type="l")
	matplot(ros[1:N1],cbind(der3.f[1:N1],(der2.f[1:N1]-der2.f[2:N])/h),type="l")
}
if(job==0.3) #check der1,2,3
{
	ros=lnf=der1.f=der2.f=der3.f=seq(from=-.95,to=.95,length=N)
	N=length(ros)
	par(mfrow=c(4,1),mar=c(3,3,1,1))
	all=matrix(nrow=N,ncol=4)
	for(i in 1:N)
	{
		all[i,1]=log(dens.ro(r,ros[i],n))
		all[i,2]=der1(r=r,ro=ros[i],n=n)
		all[i,3]=der1(r=r,ro=ros[i],n=n)
		all[i,4]=der3(r=r,ro=ros[i],n=n)
	}
	for(i in 1:4)
	{
		plot(ros,all[,i],type="l")
		segments(-2,0,2,0,col=2)
	}
}
if(job==0.4) #check cdf as a function of ro
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,3,1),cex.lab=1.5)
	x=seq(from=-.9,to=.9,by=.1);Lx=length(x)
	ros=seq(from=-.99,to=.99,length=N)
	cdfs=matrix(ncol=Lx,nrow=N)
	for(j in 1:Lx)
	for(i in 1:N)
			cdfs[i,j]=cdf.ro(r=x[j],ro=ros[i],n=n) # cdf		
	matplot(ros,cdfs,col=1,type="l",main=paste("n =",n))
}
if(job==1) #MC, MO and Olkin-Pratt estimator
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	r=r.OP=rMC=rMO=seq(from=-.999,to=.999,length=200);nr=length(r)
	for(j in 1:nr)
	{
		r.OP[j]=Re(r[j]*hypergeo(0.5,0.5,n/2-0.5,1-r[j]^2)) # Olkin - Pratt unbiased estimator
		rMO[j]=optimize(lndens.ro,r=r[j],n=n,interval=c((r[j]-1)/2,(r[j]+1)/2),maximum=T)$maximum
		rMC[j]=optimize(f=cdf.dero,r=r[j],n=n,interval=c((r[j]-1)/2,(r[j]+1)/2))$minimum		
	}
	matplot(r,cbind(r.OP-r,rMC-r,rMO-r),col=1,lwd=2,type="l",xlab="r",ylab="Difference from r")
	segments(-2,0,2,0,lty=3)		
	legend("topleft",c("Olkin-Pratt estimate","MC estimate","MO estimate"),lty=1:3,lwd=2,cex=1.2,bg="gray96")	
}
if(job==1.01) #Mode estimator
{
	ro=r
	par(mfrow=c(1,1),cex.lab=1.5,cex.main=2)	
	rMD=rep(NA,nSim)
	X=matrix(rnorm(nSim*n),ncol=n);xm=rowMeans(X)
	X2=rowSums(X^2)
	Z=matrix(rnorm(nSim*n),ncol=n)
	Y=ro*X+Z*sqrt(1-ro^2)
	Y2=rowSums(Y^2);XY=rowSums(X*Y)
	ym=rowMeans(Y)
	r=rMOD=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))
	for(i in 1:nSim)
	{
		rMD0=r[i]
		for(it in 1:10)
		{
			delta=der1(r=r[i],ro=rMD0,n=n)/der2(r=r[i],ro=rMD0,n=n)
			if(abs(delta)<0.0001) break
			rMD0=rMD0-delta			
		}
		rMOD[i]=rMD0			
	}		
	dMOD=density(rMOD)
	plot(dMOD,type="l",ylim=c(0,2.6))	
	segments(ro,-1,ro,1000,col=3)
	der=density(r)
	lines(der$x,der$y,col=2)
}
if(job==1.11) #MC estimator illustration
{
	par(mfrow=c(1,1),mar=c(5,5,1,1),cex.lab=1.5,cex.main=1.5)
	r.obs=cdf=round(seq(from=-.9,to=.9,by=.3),1);nr.obs=length(r.obs)
	N=500
	r=seq(from=-.999,to=.999,length=N)
	n=10
	pr=.3
	plot(c(-1,1),c(0,1),type="n",xlab="",ylab="",main="")
	r45x=seq(from=-1,to=1,length=N);r45y=seq(from=0,to=1,length=N)
	for(i in 1:nr.obs)
	{
		for(j in 1:N) cdf[j]=integrate(dens.ro,ro=r[j],n=n,lower=-1,upper=r.obs[i])$value
		lines(r,cdf,lwd=2)
		del=min(1-r.obs[i],r.obs[i]+1)
		ro.u=uniroot(num.int,r=r.obs[i],n=n,lower=-1*pr+(1-pr)*r.obs[i],upper=1*pr+(1-pr)*r.obs[i])$root
		cdf.i=integrate(dens.ro,ro=ro.u,n=n,lower=-1,upper=r.obs[i])$value
		points(ro.u,cdf.i,cex=1.5)	
		dd=(r45x-r)^2+(r45y-cdf)^2
		xt=r45x[dd==min(dd)];yt=cdf[dd==min(dd)]
		text(xt+.04,yt,paste("r=",r.obs[i],sep=""),adj=0,cex=1.25)
	}
}
if(job==1.12)
{	
	par(mfrow=c(1,1),mar=c(5.5,4.5,1,1),cex.lab=1.5,cex.main=1.5)
	pr=.3
	r.obs=cdf.i=ro.u=seq(from=-.9,to=.9,by=.05);nr.obs=length(r.obs)
	for(i in 1:nr.obs) 
	{
		ro.u[i]=uniroot(num.int,r=r.obs[i],n=n,lower=-1*pr+(1-pr)*r.obs[i],upper=1*pr+(1-pr)*r.obs[i])$root
		cdf.i[i]=integrate(dens.ro,ro=ro.u[i],n=n,lower=-1,upper=r.obs[i])$value
	}
	#lines(ro.u,cdf.i)
	lambda=seq(from=0.001,to=.95,length=N)
	low=up=matrix(nrow=N,ncol=2)
	r=c(-.1,.4);n=c(10,100)
	rMC=rMC.opt=rep(NA,2)
	H=matrix(ncol=2,nrow=2)
	for(i in 1:2)
	{
		for(j in 1:N)
		{
			Z1a=qnorm((1+lambda[j])/2)
			ex=exp(-2*Z1a/sqrt(n[i]-3))
			rr=(1+r[i])/(1-r[i])
			low.ro=(rr*ex-1)/(rr*ex+1)
			ex=exp(2*Z1a/sqrt(n[i]-3))
			up.ro=(rr*ex-1)/(rr*ex+1)
			qU=up.ro;qL=low.ro
			for(it in 1:10) #Optimized
			{
				itc=it
				rh1=cdf.ro(r[i],qL,n[i])-cdf.ro(r[i],qU,n[i])-lambda[j]
				rh2=cdf.dero(r[i],qL,n[i])-cdf.dero(r[i],qU,n[i])
				H[1,1]=cdf.dero(r[i],qL,n[i]);H[1,2]=-cdf.dero(r[i],qU,n[i])
				H[2,1]=num.int(r[i],qL,n[i]);H[2,2]=-num.int(r[i],qU,n[i])
				delta=solve(H)%*%c(rh1,rh2)
				q12=c(qL,qU)-delta			
				if(max(abs(delta))<10^-4) break
				qL=q12[1];qU=q12[2]			
			}
			low[j,i]=qL;up[j,i]=qU
		}
#		rMC[i]=uniroot(num.int,r=r[i],n=n[i],lower=-1*pr+(1-pr)*r[i],upper=1*pr+(1-pr)*r[i])$root
		rMC[i]=uniroot(num.int,r=r[i],n=n[i],lower=-.999,upper=.999)$root
		rMC.opt[i]=optimize(cdf.dero,r=r[i],n=n[i],lower=-.999,upper=.999)$minimum # cdf
	}
	print(rMC)
	matplot(lambda,cbind(low,up),type="l",lty=c(1,2,1,2),lwd=2,col=1,xlab="Confidence level, lambda",ylab="",main="")
	mtext(side=2,"Shortest CI",line=3,cex=1.4)
	segments(rep(-2,2),rMC,rep(.6,2),rMC,lty=1:2);text(.61,rMC,paste("MC estimate =",round(rMC,3),sep=""),adj=0,cex=1.25)	
	points(c(0,0),rMC,cex=2.5);points(c(0,0),rMC.opt,pch=3,cex=2.5)
	legend("bottomleft",c(paste("r =",r[1],", n =",n[1]),paste("r =",r[2],", n =",n[2])),lty=1:2,lwd=2,cex=1.5,bg="gray94")
}
if(job==1.2)
{
	Mr=function(ro,r,n) -(n-4)*r/(1-r^2)+(2*n-3)/2*ro/(1-ro*r)+ro/4/(2*n-1)*Re(hypergeo(0.5+1,0.5+1,n-0.5+1,r*ro/2+0.5))/Re(hypergeo(0.5,0.5,n-0.5,r*ro/2+0.5))
	
	r=ro.u=rM=ro.N=seq(from=-.9,to=.9,by=.1)
	nr=length(r)
	pr=.3
	for(i in 1:nr)
	{
		rM[i]=uniroot(Mr,r=r[i],n=n,lower=-.9,upper=.9)$root
		ro.u[i]=uniroot(num.int,r=r[i],n=n,lower=-1*pr+(1-pr)*r[i],upper=1*pr+(1-pr)*r[i])$root	
		
		ro.N[i]=r[i]
		for(it in 1:10)
		{
			ro.New=ro.N[i]-num.int(r=r[i],ro=ro.N[i],n=n)/den.int(r=r[i],ro=ro.N[i],n)
			if(abs(ro.New-ro.N[i])<0.0001) break
			ro.N[i]=ro.New			
		}	
		
	}
print(cbind(ro.u,ro.N))	
matplot(r,cbind(ro.u,ro.N))

}
if(job==2) #simulations
{
	nSim=100000
	n=15;ro=.7
	bw=.025
	X=matrix(rnorm(nSim*n),ncol=n);xm=rowMeans(X)
	X2=rowSums(X^2)
	Z=matrix(rnorm(nSim*n),ncol=n)
	Y=ro*X+Z*sqrt(1-ro^2)
	Y2=rowSums(Y^2);XY=rowSums(X*Y)
	ym=rowMeans(Y)
	r=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))
	cov.r=mean((ro<r+delta) & (ro>r-delta))
	
	Gr=Re(r*hypergeo(0.5,0.5,n/2-0.5,1-r^2))
	cov.Gr=mean((ro<Gr+delta) & (ro>Gr-delta))
		
	rM=r
	for(it in 1:10)
	{
		f0=Re(hypergeo(0.5,0.5,n-0.5,r*rM/2+0.5))
		f1=Re(hypergeo(0.5+1,0.5+1,n-0.5+1,r*rM/2+0.5))
		f2=Re(hypergeo(0.5+2,0.5+2,n-0.5+2,r*rM/2+0.5))
		hr=-rM*(n-1)/(1-rM^2)+r*(2*n-3)/2/(1-rM*r)+r/4/(2*n-1)*f1/f0
		t1=-(n-1)*(1+rM^2)/(1-rM^2)^2
		t2=r*(2*n-3)/2*r/(1-rM*r)^2
		derh=t1+t2+r/4/(2*n-1)*(r*9/(2*n+1)*f2/f0-r/4/(2*n+1)*f1^2/f0^2)
		rM.new=rM-hr/derh
		if(max(abs(rM-rM.new))<10^-4) break
		rM=rM.new
	}
	cov.rM=mean((ro<rM+delta) & (ro>rM-delta))
	print(c(cov.r,cov.Gr,cov.rM))
	denr=density(r,from=0,to=1,bw=bw);denGr=density(Gr,from=0,to=1,bw=bw);denrM=density(rM,from=-0,to=1,bw=bw)
	matplot(cbind(denr$x,denGr$x,denrM$x),cbind(denr$y,denGr$y,denrM$y),lwd=2,lty=1,type="l")
	segments(ro,-1,ro,1000)
	#rsim_7=moder(job=2)
	print(Sys.time())
	return(cbind(r,Gr,rM))
}
if(job==2.1)
{
	ro=0.7
	bw=0.05
	par(mfrow=c(1,2),cex.main=2)
	denr=density(rsim_7[,1],from=0,to=1,bw=bw);denGr=density(rsim_7[,2],from=0,to=1,bw=bw);denrM=density(rsim_7[,3],from=-0,to=1,bw=bw)
	matplot(cbind(denr$x,denGr$x,denrM$x),cbind(denr$y,denGr$y,denrM$y),lwd=2,lty=1:3,col=1,type="l",main="n = 7",xlab="",ylab="")
	mtext(side=1,"Density argument, r",line=3.5,cex=1.5)
	mtext(side=2,"Density",line=2.75,cex=1.5)
	segments(ro,-1,ro,1000)
	legend(0,2.55,c("Regular ","Unbiased (Olkin-Pratt)","Mode-estimator"),lwd=2,col=1,lty=1:3,cex=1.4,bg="gray96")
	
	denr=density(rsim_15[,1],from=0,to=1,bw=bw);denGr=density(rsim_15[,2],from=0,to=1,bw=bw);denrM=density(rsim_15[,3],from=-0,to=1,bw=bw)
	matplot(cbind(denr$x,denGr$x,denrM$x),cbind(denr$y,denGr$y,denrM$y),lwd=2,lty=1:3,col=1,type="l",main="n = 15",xlab="",ylab="")
	segments(ro,-1,ro,1000)
	mtext(side=1,"Density argument, r",line=3.5,cex=1.5)
	legend(0,3.1,c("Regular ","Unbiased (Olkin-Pratt)","Mode-estimator"),lwd=2,col=1,lty=1:3,cex=1.4,bg="gray96")
}
if(job==2.2)
{
	par(mfrow=c(1,2),cex.main=2)
	r=seq(from=-.9,to=.9,length=N)
	ro=.7
	for(n in c(7,15))
	{
		Gr=Re(r*hypergeo(0.5,0.5,n/2-0.5,1-r^2))
		rM=r
		for(it in 1:10)
		{
			f0=Re(hypergeo(0.5,0.5,n-0.5,r*rM/2+0.5))
			f1=Re(hypergeo(0.5+1,0.5+1,n-0.5+1,r*rM/2+0.5))
			f2=Re(hypergeo(0.5+2,0.5+2,n-0.5+2,r*rM/2+0.5))
			hr=-rM*(n-1)/(1-rM^2)+r*(2*n-3)/2/(1-rM*r)+r/4/(2*n-1)*f1/f0
			t1=-(n-1)*(1+rM^2)/(1-rM^2)^2
			t2=r*(2*n-3)/2*r/(1-rM*r)^2
			derh=t1+t2+r/4/(2*n-1)*(r*9/(2*n+1)*f2/f0-r/4/(2*n+1)*f1^2/f0^2)
			rM.new=rM-hr/derh
			if(max(abs(rM-rM.new))<10^-4) break
			rM=rM.new
		}
		matplot(r,cbind(Gr-r,rM-r),type="l",ylim=c(-0.04,0.04),col=1,lwd=3,lty=2:3,xlab="",ylab="",main=paste("n = ",n,", ro = ",ro,sep=""))
		segments(-2,0,2,0)
		mtext(side=1,"C.c. r",line=3.5,cex=1.5)
		mtext(side=2,"Difference from r",line=2.75,cex=1.5)
	}
	legend(-.8,-0.027,c("Unbiased (Olkin-Pratt)","Mode-estimator"),lwd=2,col=1,lty=2:3,cex=1.4,bg="gray96")
}
if(job==2.3) #test if cdf decreases with ro. Would be great to prove it
{
	ros=cdfrr=seq(from=-0.9,to=0.9,length=N)
	for(i in 1:N) cdfrr[i]=cdf.ro(r=r,ro=ros[i],n=n)
	plot(ros,cdfrr,type="l")
}
if(job==2.4) #check the formula for F derivative
{
	a=0.5;b=0.9;cc=2.4;u0=.5;u1=.9
	derF=function(a,b,cc,z) Re(a*b/cc*hypergeo(a+1,b+1,cc+1,z))
	di=integrate(f=derF,a=a,b=b,cc=cc,lower=u0,upper=u1)$value
	diF=Re(hypergeo(a,b,cc,z=u1)-hypergeo(a,b,cc,z=u0))
	print(c(di,diF))
}
if(job==3) #Z-Fisher, equal-tail, unbiased test 
{
	par(mfrow=c(1,1),mar=c(4,4.5,1,1))
	n=10;alpha=0.05;ro.true=0.8
	q1ZF=qnorm(alpha/2,mean=0.5*log((1+ro.true)/(1-ro.true)),sd=1/sqrt(n-3)) 
	q2ZF=qnorm(1-alpha/2,mean=0.5*log((1+ro.true)/(1-ro.true)),sd=1/sqrt(n-3))
	q1=q1Z=(exp(2*q1ZF)-1)/(exp(2*q1ZF)+1);q2=q2Z=(exp(2*q2ZF)-1)/(exp(2*q2ZF)+1) #Fisher-Z q1 and q2 
	
	q1UN=q1;q2UN=q2
	#Equal-tailed qLET and qUET
	qLET=q1;qUET=q2
	for(it in 1:10)
	{
		delta1=(cdf.ro(qLET,ro=ro.true,n=n)-alpha/2)/dens.ro(qLET,ro=ro.true,n=n)
		delta2=(cdf.ro(qUET,ro=ro.true,n=n)-(1-alpha/2))/dens.ro(qUET,ro=ro.true,n=n)
		qLET=qLET-delta1;qUET=qUET-delta2		
	}
	
	#Unbiased test q1=qL and q2=qU starting from Fisher-Z
	H=matrix(ncol=2,nrow=2)
	for(it in 1:maxit)
	{
		rh1=cdf.ro(q2,ro=ro.true,n=n)-cdf.ro(q1,ro=ro.true,n=n)-(1-alpha)
		rh2=cdf.dero(q1,ro=ro.true,n=n)-cdf.dero(q2,ro=ro.true,n=n)
		H[1,1]=-dens.ro(q1,ro=ro.true,n=n);H[1,2]=dens.ro(q2,ro=ro.true,n=n)
		H[2,1]=densder.ro(q1,ro=ro.true,n=n);H[2,2]=-densder.ro(q2,ro=ro.true,n=n)
		delta=solve(H)%*%c(rh1,rh2)
		q12=c(q1,q2)-delta
		if(max(abs(delta))<10^-7) break
		q1=q12[1];q2=q12[2]	
	}
	
	
		
	ros=powET=powUN=powZ=seq(from=0.25,to=0.99,length=N)
	Z1a=qnorm(1-alpha/2)
	roL1=(1+ro.true)-(1-ro.true)*exp(2*Z1a/sqrt(n-3))
	roL2=(1+ro.true)+(1-ro.true)*exp(2*Z1a/sqrt(n-3))
	roU1=(1+ro.true)-(1-ro.true)*exp(-2*Z1a/sqrt(n-3))
	roU2=(1+ro.true)+(1-ro.true)*exp(-2*Z1a/sqrt(n-3))
	#Compute 4 powers
	for(i in 1:N)
	{
		powUN[i]=1+cdf.ro(q1UN,ro=ros[i],n=n)-cdf.ro(q2UN,ro=ros[i],n=n)
		powET[i]=1+cdf.ro(qLET,ro=ros[i],n=n)-cdf.ro(qUET,ro=ros[i],n=n)
		powZ[i]=1+cdf.ro(roL1/roL2,ros[i],n=n)-cdf.ro(roU1/roU2,ros[i],n=n)
	}
	
	matplot(ros,cbind(powUN,powET,powZ),type="l",lty=c(1,2,3,4),ylim=c(0,.4),col=1,lwd=2,xlab="",ylab="Power")
	segments(-2,alpha,2,alpha,lty=2)
	segments(ro.true,-2,ro.true,2,lty=2)
	legend("bottomleft",c("1=Unbiased","2=Equal-tail","3=Fisher Z","4=Density level"),lty=c(1,2,3,4),pch=1,lwd=2,cex=1.75,bg="gray96")
	
#Simulations
	roSIM=.9
	X=matrix(rnorm(nSim*n),ncol=n)
	X2=rowSums(X^2);xm=rowMeans(X)
	Z=matrix(rnorm(nSim*n),ncol=n)
	Y=roSIM*X+Z*sqrt(1-roSIM^2)
	Y2=rowSums(Y^2);XY=rowSums(X*Y)
	ym=rowMeans(Y)
	r=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))		
	powUN.sim=1-mean(r>q1UN & r<q2UN);points(roSIM,powUN.sim,cex=1.75)
	powET.sim=1-mean(r>qLET & r<qUET);points(roSIM,powET.sim,cex=1.75)
	powZ.sim=1-mean(r>q1Z & r<q2Z);points(roSIM,powZ.sim,cex=1.75)	
	
}
if(job==3.1)#check p-value and the power
{
	n=10;alpha=0.05;ro.true=-0.2;ro.alt=0.1
	#simulations
	X=matrix(rnorm(nSim*n),ncol=n)
	X2=rowSums(X^2);xm=rowMeans(X)
	Z=matrix(rnorm(nSim*n),ncol=n)
	Y=ro.alt*X+Z*sqrt(1-ro.alt^2)
	Y2=rowSums(Y^2);XY=rowSums(X*Y)
	ym=rowMeans(Y)
	r=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))		
	FS=rep(0,nSim)
	for(i in 1:nSim) FS[i]=cdf.ro(r[i],ro=ro.true,n=n)
	
	
	q1ZF=qnorm(alpha/2,mean=0.5*log((1+ro.true)/(1-ro.true)),sd=1/sqrt(n-3)) 
	q2ZF=qnorm(1-alpha/2,mean=0.5*log((1+ro.true)/(1-ro.true)),sd=1/sqrt(n-3))
	q1=(exp(2*q1ZF)-1)/(exp(2*q1ZF)+1);q2=(exp(2*q2ZF)-1)/(exp(2*q2ZF)+1) #Fisher-Z q1 and q2 
	
	#Unbiased test q1=qL and q2=qU starting from Fisher-Z
	H=matrix(ncol=2,nrow=2)
	for(it in 1:maxit)
	{
		rh1=cdf.ro(q2,ro=ro.true,n=n)-cdf.ro(q1,ro=ro.true,n=n)-(1-alpha)
		rh2=cdf.dero(q1,ro=ro.true,n=n)-cdf.dero(q2,ro=ro.true,n=n)
		#print(c(it,q1,q2,rh1,rh2))
		H[1,1]=-dens.ro(q1,ro=ro.true,n=n);H[1,2]=dens.ro(q2,ro=ro.true,n=n)
		H[2,1]=densder.ro(q1,ro=ro.true,n=n);H[2,2]=-densder.ro(q2,ro=ro.true,n=n)
		delta=solve(H)%*%c(rh1,rh2)
		q12=c(q1,q2)-delta
		if(max(abs(delta))<10^-7) break
		q1=q12[1];q2=q12[2]	
		
	}
	q1UN=q1;q2UN=q2
	pv1=alpha/cdf.ro(q1UN,ro=ro.true,n=n)*FS
	pv2=alpha/(alpha-cdf.ro(q1UN,ro=ro.true,n=n))*(1-FS)
	a=cdf.ro(q1UN,ro=ro.true,n=n)/alpha
	
	pv=c(pv1[FS<a],pv2[FS>=a])
	powUNB=1-cdf.ro(q2UN,ro=ro.alt,n=n)+cdf.ro(q1UN,ro=ro.alt,n=n)
	out=as.data.frame(matrix(ncol=2,nrow=3))
	row.names(out)=c("Unbiased","Equal-tail","DL")
	names(out)=c("p-value","Power")
	out[1,]=c(mean(pv<alpha),powUNB)	
	
	
	#Equal-tailed qLET and qUET
	qLET=q1;qUET=q2
	for(it in 1:10)
	{
		delta1=(cdf.ro(qLET,ro=ro.true,n=n)-alpha/2)/dens.ro(qLET,ro=ro.true,n=n)
		delta2=(cdf.ro(qUET,ro=ro.true,n=n)-(1-alpha/2))/dens.ro(qUET,ro=ro.true,n=n)
		qLET=qLET-delta1;qUET=qUET-delta2		
	}
	pv1=alpha/cdf.ro(qLET,ro=ro.true,n=n)*FS
	pv2=alpha/(alpha-cdf.ro(qLET,ro=ro.true,n=n))*(1-FS)
	a=cdf.ro(qLET,ro=ro.true,n=n)/alpha
	
	pv=c(pv1[FS<a],pv2[FS>=a])
	powUNB=1-cdf.ro(qUET,ro=ro.alt,n=n)+cdf.ro(qLET,ro=ro.alt,n=n)
	out[2,]=c(mean(pv<alpha),powUNB)		
	
	
	# DL test	
	for(it in 1:10)
	{
		rh1=cdf.ro(r=q2,ro=ro.true,n=n)-cdf.ro(r=q1,ro=ro.true,n=n)-(1-alpha)
		rh2=lndens.roDL(r=q1,ro=ro.true,n=n)-lndens.roDL(r=q2,ro=ro.true,n=n)
		H[1,1]=-dens.ro(r=q1,ro=ro.true,n=n);H[1,2]=dens.ro(r=q2,ro=ro.true,n=n)
		H[2,1]=lndens.der1DL(r=q1,ro=ro.true,n=n);H[2,2]=-lndens.der1DL(r=q2,ro=ro.true,n=n)
		delta=solve(H)%*%c(rh1,rh2)
		q12=c(q1,q2)-delta
		if(max(abs(delta))<10^-7) break	
		q1=q12[1];q2=q12[2]	
	}
	q1.DL=q1;q2.DL=q2
	pv1=alpha/cdf.ro(q1.DL,ro=ro.true,n=n)*FS
	pv2=alpha/(alpha-cdf.ro(q1.DL,ro=ro.true,n=n))*(1-FS)
	a=cdf.ro(q1.DL,ro=ro.true,n=n)/alpha
	
	pv=c(pv1[FS<a],pv2[FS>=a])
	powUNB=1-cdf.ro(q2.DL,ro=ro.alt,n=n)+cdf.ro(q1.DL,ro=ro.alt,n=n)
	out[3,]=c(mean(pv<alpha),powUNB)		
	
	print(out)


}
if(job==4) #CI
{

	NIT=function(r,qL,qU,n)
	{
		H=matrix(ncol=2,nrow=2)
		rh1=cdf.ro(r,qL,n)-cdf.ro(r,qU,n)-lambda
		rh2=cdf.dero(r,qL,n)-cdf.dero(r,qU,n)
		H[1,1]=cdf.dero(r,qL,n);H[1,2]=-cdf.dero(r,qU,n)
		H[2,1]=num.int(r,qL,n);H[2,2]=-num.int(r,qU,n)
		delta=solve(H)%*%c(rh1,rh2)
		return(delta)	
	}

	n=10;alpha=0.05
	lambda=1-alpha
	Z1a=qnorm((1+lambda)/2)
	ro=seq(from=-.9,to=.9,by=.5);nr=length(ro)
	covprZ=covprM=widZ=widM=covprET=widET=prM=prZ=rep(0,nr)
	par(mfrow=c(1,2),cex.lab=1.5,cex.main=1.5)
	for(i in 1:nr)
	{
		print(c(Sys.time(),i))
		X=matrix(rnorm(nSim*n),ncol=n)
		X2=rowSums(X^2);xm=rowMeans(X)
		Z=matrix(rnorm(nSim*n),ncol=n)
		Y=ro[i]*X+Z*sqrt(1-ro[i]^2)
		Y2=rowSums(Y^2);XY=rowSums(X*Y)
		ym=rowMeans(Y)
		r=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))		
		ex=exp(-2*Z1a/sqrt(n-3))
		rr=(1+r)/(1-r)
		low.ro=(rr*ex-1)/(rr*ex+1)
		ex=exp(2*Z1a/sqrt(n-3))
		up.ro=(rr*ex-1)/(rr*ex+1)
		covprZ[i]=mean(ro[i]<up.ro & ro[i]>low.ro)	
		widZ[i]=mean(up.ro-low.ro)
		
		H=matrix(ncol=2,nrow=2)
		qU=up.ro;qL=low.ro
		for(isim in 1:nSim)
		{
			for(it in 1:10) #Optimized
			{
				itc=it
				#rh1=cdf.ro(r[isim],qL[isim],n)-cdf.ro(r[isim],qU[isim],n)-lambda
				#rh2=cdf.dero(r[isim],qL[isim],n)-cdf.dero(r[isim],qU[isim],n)
				#H[1,1]=cdf.dero(r[isim],qL[isim],n);H[1,2]=-cdf.dero(r[isim],qU[isim],n)
				#H[2,1]=num.int(r[isim],qL[isim],n);H[2,2]=-num.int(r[isim],qU[isim],n)
				#delta=solve(H)%*%c(rh1,rh2)
				
				delta=as.vector(try(NIT(r=r[isim],qL=qL[isim],qU=qU[isim],n),silent=T))
				if(!is.na(delta[1]))
				{
					q12=c(qL[isim],qU[isim])-delta			
					if(max(abs(delta))<10^-8) break				
					qL[isim]=q12[1];qU[isim]=q12[2]			
				}				
				else
				{
					qL[isim]=qU[isim]=NA
					break
				}
			}
			#print(c(isim,itc,qL[isim],qU[isim]))			
		}
		covprM[i]=mean(ro[i]<qU & ro[i]>qL,na.rm=T)		
		widM[i]=mean(qU-qL,na.rm=T)				
	}
	matplot(ro,cbind(covprM,covprZ),type="l",lwd=3,xlab="True correlation coefficient, ro",ylab="Coverage of ro")
	segments(-1,lambda,1,lambda,lwd=3,col=3)
	title(paste("Sample size n=",n,", nominal significance level =",lambda,", number of simulations =",nSim)) 
	
	matplot(ro,cbind(widM,widZ),type="l",lwd=3,xlab="True correlation coefficient, ro",ylab="CI width")
	title(paste("Sample size n=",n,", nominal significance level =",lambda,", number of simulations =",nSim)) 	
	
	
	#moder4.out=moder(job=4)
	print(Sys.time()-t0)
	return(list(n,nSim,(cbind(ro,covprZ,widZ,covprM,widM))))
}
if(job==4.1)
{
	alpha=0.05
	par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),cex.lab=1.5,cex.main=2)
	n=moder4.out[[1]]
	nSim=moder4.out[[2]]
	ma=moder4.out[[3]]
	ro=ma[,1]
	matplot(ro,cbind(ma[,2],ma[,4]),xlim=c(-1,1),ylim=c(.9,1),type="l",col=1,lwd=2,xlab="",ylab="",main="Coverage probability")
	segments(-2,1-alpha,2,1-alpha)
	mtext(side=1,"r",cex=2,line=3,font=5)
	mtext(side=2,"Proportion",line=3,cex=2)
	legend("topleft",c("Fisher Z-transformation","Optimized CI"),lty=1:2,lwd=2,cex=1.75,bg="gray96")
	matplot(ro,cbind(ma[,3],ma[,5]),type="l",xlim=c(-1,1),col=1,lwd=2,ylab="",main="Length of CI",xlab="")
	mtext(side=2,"Length",cex=2,line=3)
	mtext(side=1,"r",cex=2,line=3,font=5)
	legend("bottomright",c("Fisher Z-transformation","Optimized CI"),lty=1:2,lwd=2,cex=1.75,bg="gray96")

}
if(job==4.2) #4 CIs tested by simulations for coverage and width
{
	n=10;alpha=0.05
	lambda=1-alpha
	Z1a=qnorm((1+lambda)/2)
	ro=seq(from=-.9,to=.9,by=.1);nr=length(ro)
	covprZ=covpr.unb=widZ=wid.unb=covprET=widET=pr.short=wid.short=rep(0,nr)
	qLET=qUET=rep(NA,nSim)
	for(i in 1:nr)
	{
		X=matrix(rnorm(nSim*n),ncol=n)
		X2=rowSums(X^2);xm=rowMeans(X)
		Z=matrix(rnorm(nSim*n),ncol=n)
		Y=ro[i]*X+Z*sqrt(1-ro[i]^2)
		Y2=rowSums(Y^2);XY=rowSums(X*Y)
		ym=rowMeans(Y)
		r=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))	
		
		ex=exp(-2*Z1a/sqrt(n-3))
		rr=(1+r)/(1-r)
		low.ro=(rr*ex-1)/(rr*ex+1)
		ex=exp(2*Z1a/sqrt(n-3))
		up.ro=(rr*ex-1)/(rr*ex+1)
		#Fisher Z CI
		covprZ[i]=mean(ro[i]<up.ro & ro[i]>low.ro)	
		widZ[i]=mean(up.ro-low.ro)
				
		H=matrix(ncol=2,nrow=2)
		qU=qUsh=up.ro;qL=qLsh=low.ro
		for(isim in 1:nSim)
		{
		
			#Equal-tail CI
			qLET[isim]=low.ro[isim];qUET[isim]=up.ro[isim]
			for(it in 1:10)
			{
				delta1=(cdf.ro(r[isim],ro=qUET[isim],n=n)-alpha/2)/integrate(dlnd,ro=qUET[isim],n=n,lower=-1,upper=r[isim])$value
				delta2=(cdf.ro(r[isim],qLET[isim],n=n)-(1-alpha/2))/integrate(dlnd,ro=qLET[isim],n=n,lower=-1,upper=r[isim])$value
				qUET[isim]=qUET[isim]-delta1;qLET[isim]=qLET[isim]-delta2		
			}		
				
			for(it in 1:10) # Unbiased CI
			{
				rh1=cdf.ro(r[isim],qL[isim],n)-cdf.ro(r[isim],qU[isim],n)-lambda
				#rh2=cdf.dero(r[isim],qL[isim],n)-cdf.dero(r[isim],qU[isim],n)
				rh2=lndens.ro(qL[isim],r[isim],n)-lndens.ro(qU[isim],r[isim],n)
				H[1,1]=cdf.dero(r[isim],qL[isim],n);H[1,2]=-cdf.dero(r[isim],qU[isim],n)
				H[2,1]=der1(r[isim],qL[isim],n);H[2,2]=-der1(r[isim],qU[isim],n)
				delta=solve(H)%*%c(rh1,rh2)
				q12=c(qL[isim],qU[isim])-delta			
				if(max(abs(delta))<10^-4) break
#				print(c(isim,it,qL[isim],qU[isim],rh1,rh2))
				qL[isim]=q12[1];qU[isim]=q12[2]			
			}
			
			#Shortest CI
			qLET[isim]=qL[isim];qUsh[isim]=qUET[isim]
			for(it in 1:30) 
			{
				rh1=cdf.ro(r[isim],qLsh[isim],n)-cdf.ro(r[isim],qUsh[isim],n)-lambda
				rh2=cdf.dero(r[isim],qLsh[isim],n)-cdf.dero(r[isim],qUsh[isim],n)
				#print(c(it,isim,r[isim],rh1,rh2,qLsh[isim],qUsh[isim]))
				H[1,1]=cdf.dero(r[isim],qLsh[isim],n);H[1,2]=-cdf.dero(r[isim],qUsh[isim],n)
				H[2,1]=num.int(r[isim],qLsh[isim],n);H[2,2]=-num.int(r[isim],qUsh[isim],n)
				delta=solve(H)%*%c(rh1,rh2)
				q12=c(qLsh[isim],qUsh[isim])-delta			
				if(max(abs(delta))<10^-8) break
				qLsh[isim]=q12[1];qUsh[isim]=q12[2]			
			}
		}			
			
		covprZ[i]=mean(ro[i]<up.ro & ro[i]>low.ro)		
		widZ[i]=mean(up.ro-low.ro)
		covprET[i]=mean(ro[i]<qUET & ro[i]>qLET)		
		widET[i]=mean(qUET-qLET)
		
		covpr.unb[i]=mean(ro[i]<qU & ro[i]>qL)		
		wid.unb[i]=mean(qU-qL)				
		
		pr.short[i]=mean(ro[i]<qUsh & ro[i]>qLsh)		
		wid.short[i]=mean(qUsh-qLsh)				
	}
	#par(mfrow=c(1,2))
	#matplot(ro,cbind(covpr.unb,covprZ,covprET,pr.short),type="l",lwd=3,xlab="True correlation coefficient, ro",ylab="Coverage of ro")
	#segments(-1,lambda,1,lambda,lwd=3,col=3)
	#title(paste("Sample size n=",n,", nominal significance level =",lambda,", number of simulations =",nSim)) 
	
	#matplot(ro,cbind(wid.unb,widZ,widET,wid.short),type="l",lwd=3,xlab="True correlation coefficient, ro",ylab="CI width")
	#title(paste("Sample size n=",n,", nominal significance level =",lambda,", number of simulations =",nSim)) 
	
	#matplot(ro,cbind(wid.unb-wid.short,widZ-wid.short,widET-wid.short),type="l",lty=1,lwd=3,xlab="True correlation coefficient, ro",ylab="CI width")
	#moder42.out=moder(job=4.2)
	print(Sys.time()-t0)
	return(list(n,nSim,(cbind(ro,covprZ,widZ,covprET,widET,covpr.unb,wid.unb,pr.short,wid.short))))
}

if(job==4.201) #Check the power vi CI
{
	n=7;alpha=0.05
	r0=.8
	lambda=1-alpha
	Z1a=qnorm((1+lambda)/2)
	ro=pow=c(-.5,0,.25,r0,.9);nr=length(ro)
	covprZ=covprM=widZ=widM=covprET=widET=prM=prZ=rep(0,nr)
#	par(mfrow=c(1,3),cex.lab=1.5,cex.main=1.5)
	covPR.Z=covPR.U=matrix(ncol=nr,nrow=nr)
	for(i in 1:nr)
	{
		X=matrix(rnorm(nSim*n),ncol=n)
		X2=rowSums(X^2);xm=rowMeans(X)
		Z=matrix(rnorm(nSim*n),ncol=n)
		Y=ro[i]*X+Z*sqrt(1-ro[i]^2)
		Y2=rowSums(Y^2);XY=rowSums(X*Y)
		ym=rowMeans(Y)
		r=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))		
		ex=exp(-2*Z1a/sqrt(n-3))
		rr=(1+r)/(1-r)
		low.ro=(rr*ex-1)/(rr*ex+1)
		ex=exp(2*Z1a/sqrt(n-3))
		up.ro=(rr*ex-1)/(rr*ex+1)		
		
		H=matrix(ncol=2,nrow=2)
		qU=up.ro;qL=low.ro
		for(isim in 1:nSim)
		{
			for(it in 1:10) #Optimized
			{
				itc=it
				rh1=cdf.ro(r[isim],qL[isim],n)-cdf.ro(r[isim],qU[isim],n)-lambda
				rh2=lndens.ro(qL[isim],r[isim],n)-lndens.ro(qU[isim],r[isim],n)
				H[1,1]=cdf.dero(r[isim],qL[isim],n);H[1,2]=-cdf.dero(r[isim],qU[isim],n)
				H[2,1]=der1(r[isim],qL[isim],n);H[2,2]=-der1(r[isim],qU[isim],n)
				delta=solve(H)%*%c(rh1,rh2)
				q12=c(qL[isim],qU[isim])-delta			
				if(max(abs(delta))<10^-4) break
#				print(c(isim,it,qL[isim],qU[isim],rh1,rh2))
				qL[isim]=q12[1];qU[isim]=q12[2]			
			}			
		}
		pow[i]=1-mean(r0>qL & r0<qU)		
	}
	print(cbind(ro,pow))
	print(date())
}



if(job==4.21) #Unbiasedness of the CI
{
	n=10;alpha=0.05
	lambda=1-alpha
	Z1a=qnorm((1+lambda)/2)
	ro=seq(from=-.9,to=.9,length=N);nr=length(ro)
	covprZ=covprM=widZ=widM=covprET=widET=prM=prZ=rep(0,nr)
#	par(mfrow=c(1,3),cex.lab=1.5,cex.main=1.5)
	covPR.Z=covPR.U=matrix(ncol=nr,nrow=nr)
	for(i in 1:nr)
	{
		X=matrix(rnorm(nSim*n),ncol=n)
		X2=rowSums(X^2);xm=rowMeans(X)
		Z=matrix(rnorm(nSim*n),ncol=n)
		Y=ro[i]*X+Z*sqrt(1-ro[i]^2)
		Y2=rowSums(Y^2);XY=rowSums(X*Y)
		ym=rowMeans(Y)
		r=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))		
		ex=exp(-2*Z1a/sqrt(n-3))
		rr=(1+r)/(1-r)
		low.ro=(rr*ex-1)/(rr*ex+1)
		ex=exp(2*Z1a/sqrt(n-3))
		up.ro=(rr*ex-1)/(rr*ex+1)
		
		
		H=matrix(ncol=2,nrow=2)
		qU=up.ro;qL=low.ro
		for(isim in 1:nSim)
		{
			for(it in 1:10) #Optimized
			{
				itc=it
				rh1=cdf.ro(r[isim],qL[isim],n)-cdf.ro(r[isim],qU[isim],n)-lambda
				rh2=lndens.ro(qL[isim],r[isim],n)-lndens.ro(qU[isim],r[isim],n)
				H[1,1]=cdf.dero(r[isim],qL[isim],n);H[1,2]=-cdf.dero(r[isim],qU[isim],n)
				H[2,1]=der1(r[isim],qL[isim],n);H[2,2]=-der1(r[isim],qU[isim],n)
				delta=solve(H)%*%c(rh1,rh2)
				q12=c(qL[isim],qU[isim])-delta			
				if(max(abs(delta))<10^-4) break
#				print(c(isim,it,qL[isim],qU[isim],rh1,rh2))
				qL[isim]=q12[1];qU[isim]=q12[2]			
			}
		}
		for(j in 1:nr) covPR.Z[i,j]=mean(ro[j]<up.ro & ro[j]>low.ro)
		for(j in 1:nr) covPR.U[i,j]=mean(ro[j]<qU & ro[j]>qL)
		
#		matplot(ro,cbind(prM,prZ),type="b",main=paste("True ro =",ro[i]))
#		segments(ro[i],-1,ro[i],1,col=2)		
	}
	print(covPR.Z)
	print(covPR.U)
	contour(ro,ro,covPR.Z)
	contour(ro,ro,covPR.U,add=T,col=2)
}
if(job==4.3)
{
	par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),cex.main=1.5,cex.lab=1.5)
	n=moder42.out[[1]]
	ro=(moder42.out[[3]])[,1]
	covprZ=(moder42.out[[3]])[,2]
	widZ=(moder42.out[[3]])[,3]
	covprET=(moder42.out[[3]])[,4]
	widET=(moder42.out[[3]])[,5]
	covpr.unb=(moder42.out[[3]])[,6]
	wid.unb=(moder42.out[[3]])[,7]
	covpr.short=(moder42.out[[3]])[,8]
	wid.short=(moder42.out[[3]])[,9]
	cx=1.3
	matplot(ro,cbind(covprZ,covprET,covpr.unb,covpr.short),type="l",col=1,lwd=2,xlab="True correlation coefficient",ylab="Probability",main="Coverage probability")
	segments(-2,.95,2,.95)
	legend("bottomright",c("Fisher","Equal-tail","Unbiased","Shortest"),lwd=3,lty=1,pch=1:4,cex=1.5,bg="gray95")
	matplot(ro,cbind(widZ,widET,wid.unb,wid.short),type="l",col=1,lwd=3,xlab="True correlation coefficient",ylab="Length",main="CI length")	

}
if(job==4.31)
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,3,1),cex.main=1.5,cex.lab=1.5)
	n=moder42.out[[1]]
	ro=(moder42.out[[3]])[,1]
	covprZ=(moder42.out[[3]])[,2]
	widZ=(moder42.out[[3]])[,3]
	covprET=(moder42.out[[3]])[,4]
	widET=(moder42.out[[3]])[,5]
	covpr.unb=(moder42.out[[3]])[,6]
	wid.unb=(moder42.out[[3]])[,7]
	covpr.short=(moder42.out[[3]])[,8]
	wid.short=(moder42.out[[3]])[,9]
	cx=1.3
	matplot(ro,cbind(widZ-wid.short,widET-wid.short,wid.unb-wid.short),type="l",col=1,lty=1,lwd=2,ylim=c(0,.08),xlab="True correlation coefficient",ylab="Length",main="CI length")	
	points(ro,widZ-wid.short,pch=1,cex=1.5)
	points(ro,widET-wid.short,pch=2,cex=1.5)
	points(ro,wid.unb-wid.short,pch=3,cex=1.5)
	legend("bottomright",c("Fisher","Equal-tail","Unbiased"),lwd=3,lty=1,pch=1:3,cex=1.5,bg="gray95")
}
if(job==4.4) #Short CI: why covprob < .95?
{
	ro=.2
	lambda=1-alpha
	Z1a=qnorm((1+lambda)/2)
#r generation		
	X=matrix(rnorm(nSim*n),ncol=n)
	X2=rowSums(X^2);xm=rowMeans(X)
	Z=matrix(rnorm(nSim*n),ncol=n)
	Y=ro*X+Z*sqrt(1-ro^2)
	Y2=rowSums(Y^2);XY=rowSums(X*Y)
	ym=rowMeans(Y)
	r=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))		
	
	ex=exp(-2*Z1a/sqrt(n-3))
	rr=(1+r)/(1-r)
	low.ro=(rr*ex-1)/(rr*ex+1)
	ex=exp(2*Z1a/sqrt(n-3))
	up.ro=(rr*ex-1)/(rr*ex+1)
	covprZ=mean(ro<up.ro & ro>low.ro)	
	print(covprZ)		
		
	H=matrix(ncol=2,nrow=2)
	qL=low.ro;qU=up.ro
	for(isim in 1:nSim)
	{
		for(it in 1:100) #Optimized
		{
			rh1=cdf.ro(r[isim],qL[isim],n)-cdf.ro(r[isim],qU[isim],n)-lambda
			rh2=cdf.dero(r[isim],qL[isim],n)-cdf.dero(r[isim],qU[isim],n)				
			H[1,1]=cdf.dero(r[isim],qL[isim],n);H[1,2]=-cdf.dero(r[isim],qU[isim],n)
			H[2,1]=num.int(r[isim],qL[isim],n);H[2,2]=-num.int(r[isim],qU[isim],n)
			delta=solve(H)%*%c(rh1,rh2)
			q12=c(qL[isim],qU[isim])-delta			
			if(max(abs(delta))<10^-10) break
			qL[isim]=q12[1];qU[isim]=q12[2]			
		}		
	}
	
	pow=mean(ro>qL & ro<qU,na.rm=T)		
	print(Sys.time()-t0)
	return(pow)
}
if(job==5)
{
	ro2=r^2
	X=matrix(rnorm(nSim*n),ncol=n);xm=rowMeans(X)
	X2=rowSums(X^2)
	Z=matrix(rnorm(nSim*n),ncol=n)
	Y=r*X+Z*sqrt(1-r^2)
	Y2=rowSums(Y^2);XY=rowSums(X*Y)
	ym=rowMeans(Y)
	r=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))
	r2=r^2
	r2=sort(r2)
	plot(r2,(1:nSim)/nSim,xlim=c(0,1),ylim=c(0,1),type="s")
	x=seq(from=0,to=1,length=10000)
	sqx=sqrt(x)
	T1=0.5*log((1+sqx)/(1-sqx));T2=0.5*log((1-sqx)/(1+sqx))
	sqro=sqrt(ro2)
	R=0.5*log((1+sqro)/(1-sqro))
	Fx=pnorm(sqrt(n-3)*(T1-R))-pnorm(sqrt(n-3)*(T2-R))
	lines(x,Fx,col=2)		

}
if(job==6) #vectorized r cdf
{
	
	ro=cdf=cdf.v=seq(from=-.9,to=.2,length=N)
	r=seq(from=-.8,to=.5,length=N)
	for(i in 1:length(ro)) cdf[i]=cdf.ro(r=r[i],ro=ro[i],n)
	print(Sys.time())
	cdf.v=cdf.ro.vect(r=r,ro=ro,n)
	return(Sys.time())
	#print(cbind(ro,r,cdf,cdf.v))
	
	
	
	#ret=Vectorize(integrate(dens.ro,ro=ro,n=n,lower=-1,upper=r)$value})

}
if(job==7) #roMC, roMO and Equal-tail CI starting from Fisher Z-transformation
{
	t0=Sys.time()
	print(paste("n =",n," r =",r))
	roMC=optimize(f=cdf.dero,r=r,n=n,interval=c((r-1)/2,(r+1)/2))$minimum
	roMO=optimize(f=lndens.ro,r=r,n=n,interval=c((r-1)/2,(r+1)/2),maximum=T)$maximum
	print(paste("roMC =",round(roMC,5),", roMO =",round(roMO,5),sep=""))
	
	Z1a=qnorm((1+lambda)/2)
	ex=exp(-2*Z1a/sqrt(n-3))
	rr=(1+r)/(1-r)
	low.ZF=(rr*ex-1)/(rr*ex+1)
	ex=exp(2*Z1a/sqrt(n-3))
	up.ZF=(rr*ex-1)/(rr*ex+1)
	#print(c(low.ZF,up.ZF))
	et.L=low.ZF
	for(it in 1:maxit)
	{
		LHS=cdf.ro(r,et.L,n)-(1-alpha/2)
	#	print(c(it,LHS,et.L))
		delta=LHS/cdf.dero(r,et.L,n)
		if(abs(delta)<eps) break
		et.L=et.L-delta	
	}
	et.U=up.ZF
	for(it in 1:maxit)
	{
		LHS=cdf.ro(r,et.U,n)-alpha/2
	#	print(c(it,LHS,et.U))
		delta=LHS/cdf.dero(r,et.U,n)
		if(abs(delta)<eps) break
		et.U=et.U-delta	
	}
	print(paste(round(100*lambda),"% confidence intervals"))
	out=as.data.frame(cbind(c(low.ZF,up.ZF),c(et.L,et.U)),row.names=c("Lower limt","Upper limit"))
	names(out)=c("Fisher-Z","Equal-tail")
	print(out)
}
if(job==8)
{
	par(mfrow=c(1,1))
	ro=seq(from=-.9,to=.9,by=.1);nr=length(ro)
	cF=wF=cET=rep(NA,nr)
	for(i in 1:nr)
	{
		cF[i]=integrate(covprFisher_fun,ro=ro[i],n=n,lower=-1,upper=1)$value
		cET[i]=integrate(covprET_fun,ro=ro[i],n=n,lower=-1,upper=1)$value
		cET[i]=integrate(covprSHORT_fun,ro=ro[i],n=n,lower=-1,upper=1)$value		
		print(c(i,ro[i],cF[i],cET[i]))
		points(ro[i],cF[i]);points(ro[i],cET[i],pch=2)
	}
	print(cbind(cF,cET))
	plot(ro,cF,ylim=c(.9,1),type="o")	
}
if(job==9)
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	r=seq(from=-.999,to=.999,length=N);nr=length(r)
	L=U=qLET=qUET=qL=qU=rep(NA,nr)
	Z1a=qnorm(1-alpha/2)
	ex=exp(-2*Z1a/sqrt(n-3))	
	rr=(1+r)/(1-r)
	low.ro=(rr*ex-1)/(rr*ex+1)
	ex=exp(2*Z1a/sqrt(n-3))
	up.ro=(rr*ex-1)/(rr*ex+1)
	matplot(r,cbind(low.ro,up.ro),col=1,lty=1,lwd=2,type="l",xlab="Observed value, r",ylab="Lower and upper limits of the 95% CI")
	segments(-1,-1,1,1)
	text(.6,-.75,paste("n =",n),font=3,cex=2)
	legend("topleft",c("Fisher","Equal-tail","Short (MC)","Unbiased (MO)"),lty=1:4,lwd=2,cex=1.4,bg="gray96")
	for(i in 1:nr)
	{
		#Equal-tail
		qLET[i]=low.ro[i];qUET[i]=up.ro[i]
		for(it in 1:10)
		{
			delta1=(cdf.ro(r[i],ro=qUET[i],n=n)-alpha/2)/integrate(dlnd,ro=qUET[i],n=n,lower=-1,upper=r[i])$value
			delta2=(cdf.ro(r[i],qLET[i],n=n)-(1-alpha/2))/integrate(dlnd,ro=qLET[i],n=n,lower=-1,upper=r[i])$value
			qUET[i]=qUET[i]-delta1;qLET[i]=qLET[i]-delta2		
		}		
		#points(r[i],qUET[i],pch=2);points(r[i],qLET[i],pch=2)
				
		#Shortest CI
		L[i]=low.ro[i];U[i]=up.ro[i]
		H=matrix(ncol=2,nrow=2)
		for(it in 1:100) 
		{
		itc=it
			rh1=cdf.ro(r[i],L[i],n)-cdf.ro(r[i],U[i],n)-(1-alpha)
			rh2=cdf.dero(r[i],L[i],n)-cdf.dero(r[i],U[i],n)
			H[1,1]=cdf.dero(r[i],L[i],n);H[1,2]=-cdf.dero(r[i],U[i],n)
			H[2,1]=num.int(r[i],L[i],n);H[2,2]=-num.int(r[i],U[i],n)
			delta=solve(H)%*%c(rh1,rh2)
			q12=c(L[i],U[i])-delta			
			if(max(abs(delta))<10^-6) break
			L[i]=q12[1];U[i]=q12[2]			
		}		
		#print(c(3,r[i],itc,max(abs(delta))))
		#points(r[i],L[i]);points(r[i],U[i])
				
		# Unbiased CI
		qL[i]=low.ro[i];qU[i]=up.ro[i]
		for(it in 1:100)
		{
		itc=it
			rh1=cdf.ro(r[i],qL[i],n)-cdf.ro(r[i],qU[i],n)-lambda
			rh2=lndens.ro(qL[i],r[i],n)-lndens.ro(qU[i],r[i],n)
			H[1,1]=cdf.dero(r[i],qL[i],n);H[1,2]=-cdf.dero(r[i],qU[i],n)
			H[2,1]=der1(r[i],qL[i],n);H[2,2]=-der1(r[i],qU[i],n)
			delta=solve(H)%*%c(rh1,rh2)
			q12=c(qL[i],qU[i])-delta			
			if(max(abs(delta))<10^-6) break
			qL[i]=q12[1];qU[i]=q12[2]		
		}
		#points(r[i],qU[i],pch=2);points(r[i],qL[i],pch=3)	
	
	}
	lines(r,qLET,lty=2,lwd=2);lines(r,qUET,lty=2,lwd=2)
	lines(r,L,lty=3,lwd=2);lines(r,U,lty=3,lwd=2)
	lines(r,qL,lty=4,lwd=2);lines(r,qU,lty=4,lwd=2)	
}

print(Sys.time()-t0)
}
moder42.out <-
list(10, 500, structure(c(-0.90000000000000002, -0.80000000000000004, 
-0.69999999999999996, -0.59999999999999998, -0.5, -0.40000000000000002, 
-0.29999999999999999, -0.20000000000000001, -0.10000000000000001, 
0, 0.10000000000000001, 0.20000000000000001, 0.29999999999999999, 
0.40000000000000002, 0.5, 0.59999999999999998, 0.69999999999999996, 
0.80000000000000004, 0.90000000000000002, 0.96599999999999997, 
0.94599999999999995, 0.94599999999999995, 0.94199999999999995, 
0.94199999999999995, 0.95399999999999996, 0.95399999999999996, 
0.95999999999999996, 0.94399999999999995, 0.94599999999999995, 
0.95799999999999996, 0.93400000000000005, 0.94999999999999996, 
0.92600000000000005, 0.93000000000000005, 0.94999999999999996, 
0.94799999999999995, 0.93999999999999995, 0.93999999999999995, 
0.365522744315946, 0.60026234587217597, 0.768945135426925, 0.87854928621371398, 
0.98678228377069199, 1.05275888386662, 1.1059356925000201, 1.1392530670166401, 
1.1548243272865899, 1.1584468812039499, 1.1619716819006001, 1.12952552310813, 
1.1083706569884699, 1.03803544841365, 0.97730049314564105, 0.90004567697218096, 
0.77186037432791099, 0.606151726477969, 0.35305399522643, 0.95399999999999996, 
0.94599999999999995, 0.94399999999999995, 0.94199999999999995, 
0.94399999999999995, 0.95399999999999996, 0.95399999999999996, 
0.96199999999999997, 0.94599999999999995, 0.94599999999999995, 
0.95799999999999996, 0.93200000000000005, 0.95199999999999996, 
0.92000000000000004, 0.92800000000000005, 0.93799999999999994, 
0.94199999999999995, 0.93000000000000005, 0.94199999999999995, 
0.37989930812213601, 0.60646138279461104, 0.76385542378043803, 
0.86482731366296695, 0.96302463802220495, 1.0223724961476699, 
1.06966887237847, 1.09903157943476, 1.11258921399169, 1.1152298668591201, 
1.11806476619844, 1.08877575151696, 1.06972899409016, 1.00597407927198, 
0.95091192286415804, 0.88031283473556499, 0.76282481671859104, 
0.60725251952835502, 0.36401370533265198, 0.95399999999999996, 
0.94599999999999995, 0.94399999999999995, 0.93999999999999995, 
0.94199999999999995, 0.95399999999999996, 0.95399999999999996, 
0.96199999999999997, 0.94599999999999995, 0.94599999999999995, 
0.95799999999999996, 0.93200000000000005, 0.94999999999999996, 
0.92000000000000004, 0.92800000000000005, 0.93799999999999994, 
0.94199999999999995, 0.93000000000000005, 0.94199999999999995, 
0.37691222654741702, 0.60286267083880896, 0.76038089017440103, 
0.86161381780044299, 0.96030440048771004, 1.0200794047654, 1.0678931015741899, 
1.09765418560409, 1.11166516093987, 1.1146548584023399, 1.1178942361694599, 
1.0888483484032501, 1.07006588245056, 1.00644298180721, 0.95151935776187502, 
0.88096048929474602, 0.76344921513889397, 0.60776407965364498, 
0.36430037167717599, 0.96399999999999997, 0.93600000000000005, 
0.92600000000000005, 0.91800000000000004, 0.92200000000000004, 
0.92400000000000004, 0.93000000000000005, 0.93600000000000005, 
0.90400000000000003, 0.90200000000000002, 0.92000000000000004, 
0.89200000000000002, 0.92200000000000004, 0.89400000000000002, 
0.90200000000000002, 0.90800000000000003, 0.93799999999999994, 
0.92200000000000004, 0.93799999999999994, 0.33012014126870898, 
0.54853445109210897, 0.71001322841748404, 0.81620078056443601, 
0.92310635102526295, 0.98925819680533, 1.04319220405407, 1.0774055960692099, 
1.09331138703746, 1.0975024594316101, 1.1008565750372601, 1.0676325511318201, 
1.04567024712659, 0.97517239998942395, 0.914034771576542, 0.83781041686639701, 
0.71251819988676601, 0.55472362063776204, 0.318920933250461), .Dim = c(19L, 
9L), .Dimnames = list(NULL, c("ro", "covprZ", "widZ", "covprET", 
"widET", "covpr.unb", "wid.unb", "pr.short", "wid.short"))))
