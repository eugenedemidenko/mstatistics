expTCI <-
function(job=1,th.tru=0.5,n=3,N=100,alpha=0.05,maxit=10,nSim=100000)
{
dump("expTCI","c:\\Projects\\Mode\\expTCI.r")
if(job==1)
{
	tU=tUET=qgamma(1-alpha/2,shape=n)*th.tru
	tL=tLET=qgamma(alpha/2,shape=n)*th.tru
	H=matrix(0,ncol=2,nrow=2)
	for(it in 1:maxit)
	{

		rhs1=pgamma(tU/th.tru,shape=n)-pgamma(tL/th.tru,shape=n)-(1-alpha)
		rhs2=th.tru*n*log(tL/tU)-(tL-tU)
		H[1,1]=-1/th.tru*dgamma(tL/th.tru,shape=n)
		H[1,2]=1/th.tru*dgamma(tU/th.tru,shape=n)
		H[2,1]=th.tru*n/tL-1;H[2,2]=-th.tru*n/tU+1
		iH=solve(H)
		detH=H[1,1]*H[2,2]-H[1,2]*H[2,1]
		iH[1,1]=H[2,2];iH[1,2]=-H[1,2];iH[2,1]=-H[2,1];iH[2,2]=H[1,1]
		iH=iH/detH
		delta=iH%*%c(rhs1,rhs2)
		if(max(abs(delta))<0.00001) break
		tL=tL-delta[1];tU=tU-delta[2]
		#print(c(it,tL,tU,rhs1,rhs2))
	}

	par(mfrow=c(1,2))
	ths=pow=powET=seq(from=.1,to=1.5,length=N)
	for(i in 1:N)
	{
		X=matrix(rexp(n*nSim,rate=1/ths[i]),ncol=n)
		S=rowSums(X)
		pow[i]=mean(S>tU | S<tL)
		powET[i]=mean(S>tUET | S<tLET)
	}
	matplot(ths,cbind(pow,powET),col=1,type="l",ylim=c(0,1))
	segments(-1,alpha,1000,alpha,col=2)	

	ths=cov.pr=cov.prET=seq(from=.1,to=1.5,length=N)
	for(i in 1:N)
	{
		X=matrix(rexp(n*nSim,rate=1/ths[i]),ncol=n)
		S=rowSums(X)
		tL=S/qgamma(1-alpha/2,shape=n)
		tU=S/qgamma(alpha/2,shape=n)
		cov.prET[i]=mean(ths[i]<tU & ths[i]>tL)
		#for(it in 1:maxit)
		for(it in 1:1)
		{
			rhs1=pgamma(S/tL,shape=n)-pgamma(S/tU,shape=n)-(1-alpha)
			rhs2=n*th.tru*log(tL/tU)+S*(1/tL-1/tU)
			H11=-S/tL^2*dgamma(S/tL,shape=n)
			H12=S/tU^2*dgamma(S/tU,shape=n)
			H21=n/tL-S/tL^2
			H22=-n/tU+S/tU^2
			detH=H11*H22-H12*H21
			iH11=H22;iH12=-H12;iH21=-H21;iH22=H11
			delta1=(iH11*rhs1+iH12*rhs2)/detH
			delta2=(iH21*rhs1+iH22*rhs2)/detH
			if(max(abs(delta1))+max(abs(delta2))<0.00001) break
			tL=tL-delta1;tU=tU-delta2
			#print(it);print(tL);print(tU);print(rhs1);print(rhs2)
		}	
		cov.pr[i]=mean(ths[i]<tU & ths[i]>tL)
	}
	matplot(ths,cbind(cov.pr,cov.prET),type="l",ylim=c(0,1))
	segments(-1,1-alpha,1000,1-alpha,col=3)	
}
if(job==1.9)
{
	par(mfrow=c(1,1))
	X=matrix(rexp(n*nSim,rate=1/th.tru),ncol=n)
	S=sort(rowSums(X))
	plot(S,(1:nSim)/nSim,col=1,type="l",ylim=c(0,1))
	lines(S,pchisq(2*S/th.tru,df=2*n),col=2)
}
if(job==2)
{
	tU=tUET=th.tru*qchisq(1-alpha/2,df=2*n)/2
	tL=tLET=th.tru*qchisq(alpha/2,df=2*n)/2
	H=matrix(0,ncol=2,nrow=2)
	for(it in 1:maxit)
	{

		rhs1=pchisq(2*tU/th.tru,df=2*n)-pchisq(2*tL/th.tru,df=2*n)-(1-alpha)
		rhs2=n*th.tru*log(tU/tL)-(tU-tL)
		H[1,1]=-2/th.tru*dchisq(2*tL/th.tru,df=2*n)
		H[1,2]=2/th.tru*dchisq(2*tU/th.tru,df=2*n)
		H[2,1]=-th.tru*n/tL+1;H[2,2]=th.tru*n/tU-1
		iH=solve(H)
		detH=H[1,1]*H[2,2]-H[1,2]*H[2,1]
		iH[1,1]=H[2,2];iH[1,2]=-H[1,2];iH[2,1]=-H[2,1];iH[2,2]=H[1,1]
		iH=iH/detH
		delta=iH%*%c(rhs1,rhs2)
		if(max(abs(delta))<0.00001) break
		tL=tL-delta[1];tU=tU-delta[2]
		print(c(it,tL,tU,rhs1,rhs2))
	}

	par(mfrow=c(1,1))
	ths=pow=powET=seq(from=.2,to=1,length=N)
	for(i in 1:N)
	{
		X=matrix(rexp(n*nSim,rate=1/ths[i]),ncol=n)
		S=rowSums(X)
		pow[i]=mean(S>tU | S<tL)
		powET[i]=mean(S>tUET | S<tLET)
	}
	matplot(ths,cbind(powET,pow),col=1:2,lty=1,type="l")
	segments(-1,alpha,1000,alpha,col=2)	
	segments(th.tru,-1,th.tru,1,lty=3)
	
	pow.an.ET=1+pchisq(2*tLET/ths,df=2*n)-pchisq(2*tUET/ths,df=2*n)
	lines(ths,pow.an.ET,lty=1)
	pow.an=1+pchisq(2*tL/ths,df=2*n)-pchisq(2*tU/ths,df=2*n)
	lines(ths,pow.an,lty=2)

}	
if(job==3)
{
	n=2
	tU=tUET=th.tru*qchisq(1-alpha/2,df=2*n)/2
	tL=tLET=th.tru*qchisq(alpha/2,df=2*n)/2
	H=matrix(0,ncol=2,nrow=2)
	for(it in 1:maxit)
	{

		rhs1=pchisq(2*tU/th.tru,df=2*n)-pchisq(2*tL/th.tru,df=2*n)-(1-alpha)
		rhs2=n*th.tru*log(tU/tL)-(tU-tL)
		H[1,1]=-2/th.tru*dchisq(2*tL/th.tru,df=2*n)
		H[1,2]=2/th.tru*dchisq(2*tU/th.tru,df=2*n)
		H[2,1]=-th.tru*n/tL+1;H[2,2]=th.tru*n/tU-1
		iH=solve(H)
		detH=H[1,1]*H[2,2]-H[1,2]*H[2,1]
		iH[1,1]=H[2,2];iH[1,2]=-H[1,2];iH[2,1]=-H[2,1];iH[2,2]=H[1,1]
		iH=iH/detH
		delta=iH%*%c(rhs1,rhs2)
		if(max(abs(delta))<0.00001) break
		tL=tL-delta[1];tU=tU-delta[2]
		#print(c(it,tL,tU,rhs1,rhs2))
	}

	par(mfrow=c(1,1),mar=c(4,4.5,1,1),cex.lab=1.5)
	ths=seq(from=.2,to=1,length=N)	
	
	pow.an.ET=1+pchisq(2*tLET/ths,df=2*n)-pchisq(2*tUET/ths,df=2*n)
	pow.an=1+pchisq(2*tL/ths,df=2*n)-pchisq(2*tU/ths,df=2*n)
	matplot(ths,cbind(pow.an.ET,pow.an),col=1,lty=1:2,lwd=2,type="l",xlab="",ylab="Power")
	segments(-1,alpha,1000,alpha,lty=3)	
	segments(th.tru,-1,th.tru,1,lty=3)
	legend("topleft",c("Biased equal-tailed test","Unbiased unequal-tailed test"),col=1,lwd=2,lty=1:2,cex=1.5,bg="gray94")  
	text(.39,.17,paste("n =",n),font=4,cex=1.5)
	
	tL2=dchisq(2*tL/th.tru,df=2*n)*(4*tL/th.tru^3+4*tL^2/th.tru^4*((n-1)/2/tL*th.tru-tL/th.tru))
	tU2=dchisq(2*tU/th.tru,df=2*n)*(4*tU/th.tru^3+4*tU^2/th.tru^4*((n-1)/2/tU*th.tru-tU/th.tru))

}
}
