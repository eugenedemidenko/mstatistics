pareto <-
function(job=0,n=6,a.true=2,lambda=0.95,N=100,ss=3,maxiter=100,eps=10^-6,nSim=100000)
{
dump("pareto","c:\\Projects\\Mode\\pareto.r")
set.seed(ss)
t0=Sys.time()
if(job==0)#checking basic distributions
{
	X=runif(nSim)
	Y=1/(1-X)^(1/a.true) #Pareto distribution
	sY=sort(Y)
	par(mfrow=c(1,3),mar=c(4.5,4.5,3,1),cex.lab=1.5,cex.main=1.5)
	plot(sY,(1:nSim)/nSim,type="l",xlim=c(0,10),xlab=paste("sY, a=",a.true,sep=""),ylab="Density",main="Univariate Pareto cdf 1-1/x^a")
	x=seq(from=1,to=10,length=N)
	lines(x,1-1/x^a.true,type="l",col=2)
	
	X=1/(1-matrix(runif(nSim*n),ncol=n))^(1/a.true)
	S=rowSums(log(X))
	plot(density(S),type="l",xlab=paste("Sufficient statistic sumlog, a=",a.true,sep=""),ylab="Density",main="pdf of sumlog: dgamma")
	s=seq(from=min(S),to=max(S),length=N)
	lines(s,a.true^n/gamma(n)*s^(n-1)*exp(-a.true*s),col=2)	
	lines(s,dgamma(s,shape=n,rate=a.true),col=3)
	
	SS=sort(S)
	plot(SS,(1:nSim)/nSim,type="s",xlab=paste("Sufficient statistic sumlog, a=",a.true,sep=""),ylab="cdf",main="cdf of sumlog: pgamma")
	lines(s,pgamma(s,shape=n,rate=a.true),col=3)		
}
if(job==1)#Equal-tail CIs
{
	X=runif(n)
	Y=1/(1-X)^(1/a.true)
	S=sum(log(Y))	#sumlog
	aseq=seq(from=0.001,to=10,length=N)
	plot(aseq,pgamma(S,shape=n,rate=aseq),type="l")
	points(n/S,pgamma(S,shape=n,rate=n/S))
	segments(-1,(1-lambda)/2,1000,(1-lambda)/2,col=3)
	segments(-1,(1+lambda)/2,1000,(1+lambda)/2,col=3)
	aL=aU=n/S
	for(it in 1:maxiter)
	{
		num=pgamma(S,shape=n,rate=aL)-(1-lambda)/2
		den=pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL)
		deltaL=aL/n*num/den
		aL=aL-deltaL
		points(aL,pgamma(S,shape=n,rate=aL),col=3,cex=1.5)
		
		num=pgamma(S,shape=n,rate=aU)-(1+lambda)/2
		den=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU)
		deltaU=aU/n*num/den
		aU=aU-deltaU
		points(aU,pgamma(S,shape=n,rate=aU),col=3,pch=2,cex=1.5)
		if(abs(deltaL)+abs(deltaU)<eps) break
		print(c(it,aL,aU))
	}
}
if(job==1.1)#Simulations with equal-tail CIs
{
	X=1/(1-matrix(runif(nSim*n),ncol=n))^(1/a.true)
	S=rowSums(log(X))
	aL=aU=n/S
	for(it in 1:maxiter)
	{
		num=pgamma(S,shape=n,rate=aL)-(1-lambda)/2
		den=pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL)
		deltaL=aL/n*num/den
		aL=aL-deltaL
			
		num=pgamma(S,shape=n,rate=aU)-(1+lambda)/2
		den=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU)
		deltaU=aU/n*num/den
		aU=aU-deltaU
	
		if(max(abs(deltaL))+max(abs(deltaU))<eps) break
	}
	covprob=mean(a.true>aL & a.true<aU)
	print(covprob)
}
if(job==2) #unbiased CI
{
#Generation Pareto sumlog
	X=runif(n)
	Y=1/(1-X)^(1/a.true)
	S=sum(log(Y))	#sumlog
	#equal-tail CI
	aL=aU=n/S
	for(it in 1:maxiter)
	{
		num=pgamma(S,shape=n,rate=aL)-(1-lambda)/2
		den=pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL)
		deltaL=aL/n*num/den
		aL=aL-deltaL
				
		num=pgamma(S,shape=n,rate=aU)-(1+lambda)/2
		den=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU)
		deltaU=aU/n*num/den
		aU=aU-deltaU
		
		if(abs(deltaL)+abs(deltaU)<eps) break
		print(c(it,aL,aU))
	}
	
	for(it in 1:maxiter)
	{
		rhs1=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n,rate=aL)-lambda
		rhs2=n*log(aU/aL)-(aU-aL)*S
		print(c(it,aL,aU,rhs1,rhs2))
		a=-n/aL*(pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL))
		b=n/aU*(pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU))
		f=-n/aL+S;g=n/aU-S
		DELTA=a*g-b*f;DELTA_L=rhs1*g-rhs2*b;DELTA_U=rhs2*a-rhs1*f
		deltaL=DELTA_L/DELTA;deltaU=DELTA_U/DELTA
		if(abs(deltaL)+abs(deltaU)<eps) break
		aL=aL-deltaL;aU=aU-deltaU	
	}

}
if(job==2.1) #Simulations with unbiased CIs
{
#Generation Pareto sumlog
	X=1/(1-matrix(runif(nSim*n),ncol=n))^(1/a.true)
	S=rowSums(log(X))
	aL=aU=n/S
	#equal-tail CI	
	for(it in 1:maxiter)
	{
		num=pgamma(S,shape=n,rate=aL)-(1-lambda)/2
		den=pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL)
		deltaL=aL/n*num/den
		aL=aL-deltaL
		num=pgamma(S,shape=n,rate=aU)-(1+lambda)/2
		den=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU)
		deltaU=aU/n*num/den
		aU=aU-deltaU
		if(max(abs(deltaL))+max(abs(deltaU))<eps) break		
	}
	aL.eq=aL;aU.eq=aU
	covprob=mean(a.true>aL & a.true<aU)
	print(covprob)
	
	for(it in 1:maxiter)
	{
		rhs1=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n,rate=aL)-lambda
		rhs2=n*log(aU/aL)-(aU-aL)*S
		a=-n/aL*(pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL))
		b=n/aU*(pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU))
		f=-n/aL+S;g=n/aU-S
		DELTA=a*g-b*f;DELTA_L=rhs1*g-rhs2*b;DELTA_U=rhs2*a-rhs1*f
		deltaL=DELTA_L/DELTA;deltaU=DELTA_U/DELTA
		if(max(abs(deltaL))+max(abs(deltaU))<eps) break
		aL=aL-deltaL;aU=aU-deltaU	
	}
	covprob=mean(a.true>aL & a.true<aU)
	print(covprob)	
}
if(job==3) #MC short CI
{
#Generation Pareto sumlog
	X=runif(n)
	Y=1/(1-X)^(1/a.true)
	S=sum(log(Y))	#sumlog
	#equal-tail CI
	aL=aU=n/S
	for(it in 1:maxiter)
	{
		num=pgamma(S,shape=n,rate=aL)-(1-lambda)/2
		den=pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL)
		deltaL=aL/n*num/den
		aL=aL-deltaL
				
		num=pgamma(S,shape=n,rate=aU)-(1+lambda)/2
		den=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU)
		deltaU=aU/n*num/den
		aU=aU-deltaU
		
		if(abs(deltaL)+abs(deltaU)<eps) break
		print(c(it,aL,aU))
	}
	
	for(it in 1:maxiter)
	{
		rhs1=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n,rate=aL)-lambda
		rhs2=n/aU*(pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU))-n/aL*(pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL))
		print(c(it,aL,aU,rhs1,rhs2))
		a=-n/aL*(pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL))
		b=n/aU*(pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU))
		f=-n/aL^2*(n*pgamma(S,shape=n,rate=aL)-(2*n+1)*pgamma(S,shape=n+1,rate=aL)+(n+1)*pgamma(S,shape=n+2,rate=aL))
		g=n/aU^2*(n*pgamma(S,shape=n,rate=aU)-(2*n+1)*pgamma(S,shape=n+1,rate=aU)+(n+1)*pgamma(S,shape=n+2,rate=aU))
		DELTA=a*g-b*f;DELTA_L=rhs1*g-rhs2*b;DELTA_U=rhs2*a-rhs1*f
		deltaL=DELTA_L/DELTA;deltaU=DELTA_U/DELTA
		if(abs(deltaL)+abs(deltaU)<eps) break
		aL=aL-deltaL;aU=aU-deltaU	
	}
}
if(job==3.1) #Simulations with MC short CI
{
#Generation Pareto sumlog
	X=1/(1-matrix(runif(nSim*n),ncol=n))^(1/a.true)
	S=rowSums(log(X))
	#equal-tail CI
	aL=aU=n/S
	for(it in 1:maxiter)
	{
		num=pgamma(S,shape=n,rate=aL)-(1-lambda)/2
		den=pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL)
		deltaL=aL/n*num/den
		aL=aL-deltaL
				
		num=pgamma(S,shape=n,rate=aU)-(1+lambda)/2
		den=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU)
		deltaU=aU/n*num/den
		aU=aU-deltaU
		
		if(max(abs(deltaL))+max(abs(deltaU))<eps) break
		#print(c(it,aL,aU))
	}
	
	for(it in 1:maxiter)
	{
		rhs1=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n,rate=aL)-lambda
		rhs2=n/aU*(pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU))-n/aL*(pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL))
		#print(c(it,aL,aU,rhs1,rhs2))
		a=-n/aL*(pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL))
		b=n/aU*(pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU))
		f=-n/aL^2*(n*pgamma(S,shape=n,rate=aL)-(2*n+1)*pgamma(S,shape=n+1,rate=aL)+(n+1)*pgamma(S,shape=n+2,rate=aL))
		g=n/aU^2*(n*pgamma(S,shape=n,rate=aU)-(2*n+1)*pgamma(S,shape=n+1,rate=aU)+(n+1)*pgamma(S,shape=n+2,rate=aU))
		DELTA=a*g-b*f;DELTA_L=rhs1*g-rhs2*b;DELTA_U=rhs2*a-rhs1*f
		deltaL=DELTA_L/DELTA;deltaU=DELTA_U/DELTA
		if(max(abs(deltaL))+max(abs(deltaU))<eps) break
		aL=aL-deltaL;aU=aU-deltaU	
	}
	covprob=mean(a.true>aL & a.true<aU)
	print(covprob)	
}
if(job==4)# Comparison of the length of 3 CIs as a function of n via simulations
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	nes=seq(from=3,to=10,by=2);Lnes=length(nes)
	wid=matrix(ncol=3,nrow=Lnes)
	for(inn in 1:Lnes)
	{
		n=nes[inn]
		#Generation Pareto sumlog
		X=1/(1-matrix(runif(nSim*n),ncol=n))^(1/a.true)
		S=rowSums(log(X))
		#equal-tail CI
		aL=aU=n/S
		for(it in 1:maxiter)
		{
			num=pgamma(S,shape=n,rate=aL)-(1-lambda)/2
			den=pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL)
			deltaL=aL/n*num/den
			aL=aL-deltaL
				
			num=pgamma(S,shape=n,rate=aU)-(1+lambda)/2
			den=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU)
			deltaU=aU/n*num/den
			aU=aU-deltaU
		
			if(max(abs(deltaL))+max(abs(deltaU))<eps) break
			#print(c(it,aL,aU))
		}
		wid[inn,1]=mean(aU-aL)
		
		for(it in 1:maxiter)
		{
			rhs1=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n,rate=aL)-lambda
			rhs2=n*log(aU/aL)-(aU-aL)*S
			a=-n/aL*(pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL))
			b=n/aU*(pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU))
			f=-n/aL+S;g=n/aU-S
			DELTA=a*g-b*f;DELTA_L=rhs1*g-rhs2*b;DELTA_U=rhs2*a-rhs1*f
			deltaL=DELTA_L/DELTA;deltaU=DELTA_U/DELTA
			if(max(abs(deltaL))+max(abs(deltaU))<eps) break
			aL=aL-deltaL;aU=aU-deltaU	
		}
		wid[inn,2]=mean(aU-aL)
		
		for(it in 1:maxiter)
		{
			rhs1=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n,rate=aL)-lambda
			rhs2=n/aU*(pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU))-n/aL*(pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL))
			#print(c(it,aL,aU,rhs1,rhs2))
			a=-n/aL*(pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL))
			b=n/aU*(pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU))
			f=-n/aL^2*(n*pgamma(S,shape=n,rate=aL)-(2*n+1)*pgamma(S,shape=n+1,rate=aL)+(n+1)*pgamma(S,shape=n+2,rate=aL))
			g=n/aU^2*(n*pgamma(S,shape=n,rate=aU)-(2*n+1)*pgamma(S,shape=n+1,rate=aU)+(n+1)*pgamma(S,shape=n+2,rate=aU))
			DELTA=a*g-b*f;DELTA_L=rhs1*g-rhs2*b;DELTA_U=rhs2*a-rhs1*f
			deltaL=DELTA_L/DELTA;deltaU=DELTA_U/DELTA
			if(max(abs(deltaL))+max(abs(deltaU))<eps) break
			aL=aL-deltaL;aU=aU-deltaU	
		}
		wid[inn,3]=mean(aU-aL)	
	}
	matplot(nes,wid,type="o",lty=1,col=1,lwd=2,xlab="Sample size, n",ylab=paste("Length of the ",lambda*100,"% CI",sep=""))
	legend("topright",c("1: Equal-tail","2: Unbiased","3: Short"),cex=1.75,bg="grey95")
	text(4,3.5,paste("a=",a.true,"\nnSim=",nSim,sep=""),font=3,cex=1.5)
}
if(job==5)# Comparison of the length of 3 CIs as a function of a.true via simulations
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	aseq=seq(from=.5,to=2,by=.5);La=length(aseq)
	wid=matrix(ncol=3,nrow=La)
	for(inn in 1:La)
	{
		a.true=aseq[inn]
		#Generation Pareto sumlog
		X=1/(1-matrix(runif(nSim*n),ncol=n))^(1/a.true)
		S=rowSums(log(X))
		#equal-tail CI
		aL=aU=n/S
		for(it in 1:maxiter)
		{
			num=pgamma(S,shape=n,rate=aL)-(1-lambda)/2
			den=pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL)
			deltaL=aL/n*num/den
			aL=aL-deltaL
				
			num=pgamma(S,shape=n,rate=aU)-(1+lambda)/2
			den=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU)
			deltaU=aU/n*num/den
			aU=aU-deltaU
		
			if(max(abs(deltaL))+max(abs(deltaU))<eps) break
			#print(c(it,aL,aU))
		}
		wid[inn,1]=mean(aU-aL)
		
		for(it in 1:maxiter)
		{
			rhs1=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n,rate=aL)-lambda
			rhs2=n*log(aU/aL)-(aU-aL)*S
			a=-n/aL*(pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL))
			b=n/aU*(pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU))
			f=-n/aL+S;g=n/aU-S
			DELTA=a*g-b*f;DELTA_L=rhs1*g-rhs2*b;DELTA_U=rhs2*a-rhs1*f
			deltaL=DELTA_L/DELTA;deltaU=DELTA_U/DELTA
			if(max(abs(deltaL))+max(abs(deltaU))<eps) break
			aL=aL-deltaL;aU=aU-deltaU	
		}
		wid[inn,2]=mean(aU-aL)
		
		for(it in 1:maxiter)
		{
			rhs1=pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n,rate=aL)-lambda
			rhs2=n/aU*(pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU))-n/aL*(pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL))
			#print(c(it,aL,aU,rhs1,rhs2))
			a=-n/aL*(pgamma(S,shape=n,rate=aL)-pgamma(S,shape=n+1,rate=aL))
			b=n/aU*(pgamma(S,shape=n,rate=aU)-pgamma(S,shape=n+1,rate=aU))
			f=-n/aL^2*(n*pgamma(S,shape=n,rate=aL)-(2*n+1)*pgamma(S,shape=n+1,rate=aL)+(n+1)*pgamma(S,shape=n+2,rate=aL))
			g=n/aU^2*(n*pgamma(S,shape=n,rate=aU)-(2*n+1)*pgamma(S,shape=n+1,rate=aU)+(n+1)*pgamma(S,shape=n+2,rate=aU))
			DELTA=a*g-b*f;DELTA_L=rhs1*g-rhs2*b;DELTA_U=rhs2*a-rhs1*f
			deltaL=DELTA_L/DELTA;deltaU=DELTA_U/DELTA
			if(max(abs(deltaL))+max(abs(deltaU))<eps) break
			aL=aL-deltaL;aU=aU-deltaU	
		}
		wid[inn,3]=mean(aU-aL)			
	}
	matplot(aseq,wid/aseq,type="o",lty=1,col=1,lwd=2,xlab="a.true",ylab=paste("Relative length of the ",lambda*100,"% CI",sep=""))
	legend("topright",c("1: Equal-tail","2: Unbiased","3: Short"),cex=1.75,bg="grey95")
	text(1,3.5,paste("a=",a.true,"\nnSim=",nSim,sep=""),font=3,cex=1.5)
}
if(job==6)# MC estimator: Use, say, nSim=100: aMC=aMO
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	
	#Generation Pareto sumlog
	X=1/(1-matrix(runif(nSim*n),ncol=n))^(1/a.true)
	S=rowSums(log(X))
	aMO.ML=n/S
	aMC=aMO.ML
	for(it in 1:maxiter)
	{
		Gn=pgamma(S,shape=n,rate=aMC)
		Gn1=pgamma(S,shape=n+1,rate=aMC)
		Gn2=pgamma(S,shape=n+2,rate=aMC)
		Gn3=pgamma(S,shape=n+3,rate=aMC)
		num=n*Gn-(2*n+1)*Gn1+(n+1)*Gn2
		den=n^2*Gn-(3*n^2+3*n+1)*Gn1+(3*n^2+6*n+3)*Gn2-(n^2+3*n+2)*Gn3
		delta=num/den*aMC
		if(max(abs(delta))<eps) break
		aMC=aMC-delta	
	}
	plot(1:nSim,aMO.ML-aMC,ylim=c(-10^-10,10^-10),xlab="Simulation index",ylab="Difference between aMC and aMO")
}
if(job==7)# a.true=null value
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	a.alt.RR=a.true*.82;a.alt.CI=a.true*1.4
	X=1/(1-matrix(runif(nSim*n),ncol=n))^(1/a.alt.RR)
	S.RR=rowSums(log(X))
	X=1/(1-matrix(runif(nSim*n),ncol=n))^(1/a.alt.CI)
	S.CI=rowSums(log(X))
	
	a.alt=seq(from=a.true*.7,to=a.true*1.9,length=N)
	qL.et=qgamma((1-lambda)/2,shape=n,rate=a.true)
	qU.et=qgamma((1+lambda)/2,shape=n,rate=a.true)	
	pow.et=1-pgamma(qU.et,shape=n,rate=a.alt)+pgamma(qL.et,shape=n,rate=a.alt)
	plot(a.alt,pow.et,type="l",lwd=2,ylim=c(0.025,.18),xlab="Alternative Pareto parameter, a",ylab="Power")
	segments(a.alt.RR,-1,a.alt.RR,1,col="gray80");text(a.alt.RR+.05,.135,"Rejection rule",srt=90,adj=0,font=3,cex=1.25)
	segments(a.alt.CI,-1,a.alt.CI,1,col="gray80");text(a.alt.CI+.05,.135,"Confidence interval",srt=90,adj=0,font=3,cex=1.25)
	segments(-1,1-lambda,1000,1-lambda,lty=2)
	pow.RR=mean(S.RR<qL.et | S.RR>qU.et)
	points(a.alt.RR,pow.RR,cex=1.75)
	
	#equal-tail CI
	aL=aU=n/S.CI
	for(it in 1:maxiter)
	{
		num=pgamma(S.CI,shape=n,rate=aL)-(1-lambda)/2
		den=pgamma(S.CI,shape=n,rate=aL)-pgamma(S.CI,shape=n+1,rate=aL)
		deltaL=aL/n*num/den
		aL=aL-deltaL
				
		num=pgamma(S.CI,shape=n,rate=aU)-(1+lambda)/2
		den=pgamma(S.CI,shape=n,rate=aU)-pgamma(S.CI,shape=n+1,rate=aU)
		deltaU=aU/n*num/den
		aU=aU-deltaU
		
		if(max(abs(deltaL))+max(abs(deltaU))<eps) break
		#print(c(it,aL,aU))
	}
	pow.CI=mean(a.true<aL | a.true>aU)
	points(a.alt.CI,pow.CI,cex=1.75)
	
	#unbiased test
	qL.unb=qL.et;qU.unb=qU.et
	for(it in 1:maxiter)
	{
	
		rhs1=pgamma(qU.unb,shape=n,rate=a.true)-pgamma(qL.unb,shape=n,rate=a.true)-lambda
		rhs2=n/a.true*(pgamma(qU.unb,shape=n,rate=a.true)-pgamma(qU.unb,shape=n+1,rate=a.true))-n/a.true*(pgamma(qL.unb,shape=n,rate=a.true)-pgamma(qL.unb,shape=n+1,rate=a.true))
		#print(c(it,qL.unb,qU.unb,rhs1,rhs2))
		a=-n/a.true*(pgamma(qL.unb,shape=n,rate=a.true)-pgamma(qL.unb,shape=n+1,rate=a.true))
		b=n/a.true*(pgamma(qU.unb,shape=n,rate=a.true)-pgamma(qU.unb,shape=n+1,rate=a.true))
		f=-n/a.true^2*(n*pgamma(qL.unb,shape=n,rate=a.true)-(2*n+1)*pgamma(qL.unb,shape=n+1,rate=a.true)+(n+1)*pgamma(qL.unb,shape=n+2,rate=a.true))
		g=n/a.true^2*(n*pgamma(qU.unb,shape=n,rate=a.true)-(2*n+1)*pgamma(qU.unb,shape=n+1,rate=a.true)+(n+1)*pgamma(qU.unb,shape=n+2,rate=a.true))
		DELTA=a*g-b*f;DELTA_L=rhs1*g-rhs2*b;DELTA_U=rhs2*a-rhs1*f
		deltaL=DELTA_L/DELTA;deltaU=DELTA_U/DELTA
		if(max(abs(deltaL))+max(abs(deltaU))<eps) break
		qL.unb=qL.unb-deltaL;qU.unb=qU.unb-deltaU		
	}
	pow.unb=1-pgamma(qU.unb,shape=n,rate=a.alt)+pgamma(qL.unb,shape=n,rate=a.alt)
	lines(a.alt,pow.unb,lty=2,lwd=2)
	pow.RR=mean(S.RR<qL.unb | S.RR>qU.unb)
	points(a.alt.RR,pow.RR,pch=2,cex=1.75)
	
	for(it in 1:maxiter)
	{
		rhs1=pgamma(S.CI,shape=n,rate=aU)-pgamma(S.CI,shape=n,rate=aL)-lambda
		rhs2=n*log(aU/aL)-(aU-aL)*S.CI
		a=-n/aL*(pgamma(S.CI,shape=n,rate=aL)-pgamma(S.CI,shape=n+1,rate=aL))
		b=n/aU*(pgamma(S.CI,shape=n,rate=aU)-pgamma(S.CI,shape=n+1,rate=aU))
		f=-n/aL+S.CI;g=n/aU-S.CI
		DELTA=a*g-b*f;DELTA_L=rhs1*g-rhs2*b;DELTA_U=rhs2*a-rhs1*f
		deltaL=DELTA_L/DELTA;deltaU=DELTA_U/DELTA
		if(max(abs(deltaL))+max(abs(deltaU))<eps) break
		aL=aL-deltaL;aU=aU-deltaU	
	}
	pow.CI=mean(a.true<aL | a.true>aU)
	points(a.alt.CI,pow.CI,pch=2,cex=1.75)
	
	#DL test
	qL.DL=qL.et;qU.DL=qU.et
	for(it in 1:maxiter)
	{
	
		rhs1=pgamma(qU.DL,shape=n,rate=a.true)-pgamma(qL.DL,shape=n,rate=a.true)-lambda
		rhs2=(n-1)*log(qU.DL/qL.DL)-a.true*(qU.DL-qL.DL)
		#print(c(it,qL.DL,qU.DL,rhs1,rhs2))
		a=-n/a.true*(pgamma(qL.DL,shape=n,rate=a.true)-pgamma(qL.DL,shape=n+1,rate=a.true))
		b=n/a.true*(pgamma(qU.DL,shape=n,rate=a.true)-pgamma(qU.DL,shape=n+1,rate=a.true))
		f=-(n-1)/qL.DL+a.true
		g=(n-1)/qU.DL-a.true
		DELTA=a*g-b*f;DELTA_L=rhs1*g-rhs2*b;DELTA_U=rhs2*a-rhs1*f
		deltaL=DELTA_L/DELTA;deltaU=DELTA_U/DELTA
		if(max(abs(deltaL))+max(abs(deltaU))<eps) break
		qL.DL=qL.DL-deltaL;qU.DL=qU.DL-deltaU		
	}
	pow.DL=1-pgamma(qU.DL,shape=n,rate=a.alt)+pgamma(qL.DL,shape=n,rate=a.alt)
	lines(a.alt,pow.DL,lty=3,lwd=2)
	pow.RR=mean(S.RR<qL.DL | S.RR>qU.DL)
	points(a.alt.RR,pow.RR,pch=3,cex=1.75)
	legend("bottomright",c("Equal-tail","Unbiased","Density level"),col=1,lty=c(1,2,3),pch=c(1,2,3),cex=1.75,lwd=2,bg="gray95")
	
	

}
return(Sys.time()-t0)
}
