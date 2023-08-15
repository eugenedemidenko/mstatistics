unbtestr2 <-
function(n=30,p=7,ro02=c(.6,.3),r2.obs=0.4,lambda=0.95,maxit=10,N=100,ss=3,nSim=3000)
{
dump("unbtestr2","c:\\Projects\\Mode\\unbtestr2.r")
#install.packages("hypergeo")
library(hypergeo)
set.seed(ss)
t0=Sys.time()
	fr2=function(x,ro2,n,p)
	{
		fr2=gamma((n-1)/2)/gamma((n-p-1)/2)/gamma(p/2)*(1-ro2)^((n-1)/2)*(x)^((p-2)/2)*(1-x)^((n-p-3)/2)
		f0=Re(hypergeo((n-1)/2,(n-1)/2,p/2,x*ro2))
		fr2=fr2*f0 #density of r^2	
		return(fr2)
	}
	derfr2=function(x,ro2,n,p)
	{
		fr2=gamma((n-1)/2)/gamma((n-p-1)/2)/gamma(p/2)*(1-ro2)^((n-1)/2)*(x)^((p-2)/2)*(1-x)^((n-p-3)/2)
		f0=Re(hypergeo((n-1)/2,(n-1)/2,p/2,x*ro2))
		fr2=fr2*f0 #density of r^2	
		
		f1=Re(hypergeo((n+1)/2,(n+1)/2,(p+2)/2,ro2*x))
		hr=-(n-1)/2/(1-ro2)+x*(n-1)^2/2/p*f1/f0
		return(fr2*hr)		
	}

	par(mfrow=c(1,2),mar=c(4.5,4.5,4,1),cex.main=1.5,cex.lab=1.5)
	alpha=1-lambda	
	
	for(IGRAPH in 1:2)
	{
		ro2=ro02[IGRAPH]
	
	#Equal-tail
	qL.ET=ro2
	for(it in 1:maxit)
	{
		Fqr=integrate(fr2,ro2=ro2,n=n,p=p,lower=0,upper=qL.ET)$value
		delta=(Fqr-(1-lambda)/2)/fr2(x=qL.ET,ro2=ro2,n=n,p=p)
		if(abs(delta)<10^-8) break
		qL.ET=qL.ET-delta
	}
	qU.ET=(ro2+1)/2
	for(it in 1:maxit)
	{
		Fqr=integrate(fr2,ro2=ro2,n=n,p=p,lower=0,upper=qU.ET)$value
		delta=(Fqr-(1+lambda)/2)/fr2(x=qU.ET,ro2=ro2,n=n,p=p)
		if(abs(delta)<10^-8) break
		qU.ET=qU.ET-delta
	}

	qL=qL.ET;qU=qU.ET	
	H=matrix(ncol=2,nrow=2)
	for(it in 1:maxit)
	{
		eq1=integrate(fr2,ro2=ro2,n=n,p=p,lower=0,upper=qU)$value-integrate(fr2,ro2=ro2,n=n,p=p,lower=0,upper=qL)$value-(1-alpha)
		eq2=integrate(derfr2,ro2=ro2,n=n,p=p,lower=0,upper=qL)$value-integrate(derfr2,ro2=ro2,n=n,p=p,lower=0,upper=qU)$value
		#print(c(it,qL,qU,eq1,eq2))			
		H[1,1]=-fr2(qL,ro2,n,p);H[1,2]=fr2(qU,ro2,n,p)
		H[2,1]=derfr2(qL,ro2,n,p);H[2,2]=-derfr2(qU,ro2,n,p)
		delta=solve(H)%*%c(eq1,eq2)
		qL=qL-delta[1];qU=qU-delta[2]
		if(max(abs(c(eq1,eq2)))<0.0001) break		
	}		
	
#p-value
	
	Fr2=integrate(fr2,ro2=ro2,n=n,p=p,lower=0,upper=r2.obs)$value
	FqL=integrate(fr2,ro2=ro2,n=n,p=p,lower=0,upper=qL)$value
	if(Fr2<=FqL/alpha) pvalue=alpha*Fr2/FqL else pvalue=alpha*(1-Fr2)/(1-FqL)
	cat("\nObserved r2=",r2.obs,", r20=",ro2,", p-value=",round(pvalue,6),sep="")	
	cat("\nqL=",round(qL,6),", qU=",round(qU,6),sep="")	
	ro2.alt=pow=powET=powMBESS=seq(from=0,to=1,length=N)
	for(i in 1:N)
	{
		pow[i]=1+integrate(fr2,ro2=ro2.alt[i],n=n,p=p,lower=0,upper=qL)$value-integrate(fr2,ro2=ro2.alt[i],n=n,p=p,lower=0,upper=qU)$value
		powET[i]=1+integrate(fr2,ro2=ro2.alt[i],n=n,p=p,lower=0,upper=qL.ET)$value-integrate(fr2,ro2=ro2.alt[i],n=n,p=p,lower=0,upper=qU.ET)$value
	}
	#print(c(qL,qL.ET,qU,qU.ET))
	#print(cbind(ro2.alt,pow,powET))
	plot(ro2.alt,pow,type="l",xlab="Alternative ro2",ylab="Power",main=paste("n=",n,", p=",p,sep=""),lwd=3,ylim=c(0,1))
	lines(ro2.alt,powET,lty=2)
	segments(ro2,-1,ro2,1)
	text(ro2+.05,.6,paste("ro2 =",ro2),srt=90,font=3,cex=1.5)
	segments(-1,alpha,1,alpha)	
	
	#for(i in 1:N) 
	# powMBESS[i]=ss.power.R2(Population.R2=ro2.alt[i],alpha.level=alpha,p=p,Specified.N=n,Null.R2=ro2)$Actual.Power
	#lines(ro2.alt,powMBESS,lty=2)
		
	
# check with simulations
	set.seed(ss)
	A=matrix(runif((p+1)^2,min=-2,max=8),ncol=p+1)
	OM=t(A)%*%A
	vary=OM[1,1];om=OM[1,2:(p+1)]
	OMx=OM[2:(p+1),2:(p+1)]
	r2=rep(NA,nSim)
	ro2s=powsim=seq(from=0.001,to=.99,by=0.1);NROS=length(ro2s)
	for(i in 1:NROS)
	{
		OM[1,1]=t(om)%*%solve(OMx)%*%om/ro2s[i]
		mu=runif(p+1,min=-1,max=7)
		chOM=chol(OM)
		y=rep(NA,n);X=matrix(ncol=n,nrow=p)
		en=rep(1,n)
		roCI=r2=B0=rep(NA,nSim)
		for(isim in 1:nSim)
		{
			yX=matrix(rnorm(n*(p+1)),ncol=p+1,nrow=n)%*%chOM+en%*%t(mu)
			y=yX[,1]
			X=yX[,2:(p+1)]
			o=lm(y~X)
			r2[isim]=1-sum(o$residuals^2)/var(y)/(n-1)
		}
		powsim[i]=1-mean(r2>qL & r2<qU)
	}
	points(ro2s,powsim,cex=1.5,col=2,pch=16)		
}
cat("\n")
Sys.time()-t0
}
