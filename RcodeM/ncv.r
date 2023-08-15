ncv <-
function(kap0=.5,n=10,lambda=0.95,N=200,maxit=20,eps=0.001,nSim=10000,ss=6)
{
dump("ncv","c:\\Projects\\Mode\\ncv.r")
#Testing power function for an equal-tailed test
t0=Sys.time()
cdf.es=function(w,s,n,d)	2*w*pnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
cdf.der1.es=function(w,s,n,d)	-2*sqrt(n)*w*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
cdf.der2.es=function(w,s,n,d)	-2*n*w*(s*w/sqrt(n-1)-d*sqrt(n))*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
cdf.der3.es=function(w,s,n,d)	-2*n^1.5*w*((s*w/sqrt(n-1)-d*sqrt(n))^2-1)*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
f.der=function(w,s,n,d)	2*sqrt(n)*w^2/sqrt(n-1)*(s*w/sqrt(n-1)-d*sqrt(n))*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
dens=function(w,s,n,d)	2/sqrt(n-1)*w^2*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
dens.der=function(w,s,n,d)	2*sqrt(n)/sqrt(n-1)*w^2*(s*w/sqrt(n-1)-sqrt(n)*d)*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
dens.prime=function(w,s,n,d) -2/(n-1)*w^3*(s*w/sqrt(n-1)-sqrt(n)*d)*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
dens.der.prime=function(w,s,n,d) 2*sqrt(n)/(n-1)*w^3*(1+(s*w/sqrt(n-1)-sqrt(n)*d)^2)*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)

CI.ET=function(sampleES,n,lambda,maxit=100,eps=10^-7)
#the same as for effect size but returns reciprocal
{
	S=sqrt(n)*sampleES
	# Johnson and Welch (1940) approximations	
	low.es=(S-qnorm(1-alpha/2)*sqrt(1+S^2/2/n))/sqrt(n)
	up.es=(S-qnorm(alpha/2)*sqrt(1+S^2/2/n))/sqrt(n)
	for(it in 1:maxit)
	{
		#low limit
		cdf0=integrate(cdf.es,s=S,n=n,d=low.es,lower=0,upper=Inf)$value
		der1=integrate(cdf.der1.es,s=S,n=n,d=low.es,lower=0,upper=Inf)$value
		delta1=(cdf0-(1-alpha/2))/der1
		low.es=low.es-delta1
		#upper limit
		cdf0=integrate(cdf.es,s=S,n=n,d=up.es,lower=0,upper=Inf)$value
		der1=integrate(cdf.der1.es,s=S,n=n,d=up.es,lower=0,upper=Inf)$value
		delta2=(cdf0-alpha/2)/der1
		up.es=up.es-delta2
		if((abs(delta1)+abs(delta2))<eps) break
		#print(c(it,low.es,up.es))
	}
	return(c(1/up.es,1/low.es))
}
CI.unb=function(sampleES,n,lambda,low.unb,up.unb,maxit=100,eps=10^-7)
#the same as for effect size but returns reciprocals
#low.unb and up.unb must be low and up ET for CV (kappa)
{
	S=sqrt(n)*sampleES
	H=matrix(ncol=2,nrow=2)
	low.es=1/up.unb;up.es=1/low.unb
	for(it in 1:maxit)
	{
		#print(c(it,low.es,up.es))
		rhs1=pt(S,df=n-1,ncp=sqrt(n)*low.es)-pt(S,df=n-1,ncp=sqrt(n)*up.es) - (1-alpha)
		H[1,1]=integrate(cdf.der1.es,s=S,n=n,d=low.es,lower=0,upper=Inf)$value
		H[1,2]=-integrate(cdf.der1.es,s=S,n=n,d=up.es,lower=0,upper=Inf)$value
		rhs2=integrate(dens,s=S,n=n,d=low.es,lower=0,upper=Inf)$value-integrate(dens,s=S,n=n,d=up.es,lower=0,upper=Inf)$value
		H[2,1]=integrate(dens.der,s=S,n=n,d=low.es,lower=0,upper=Inf)$value
		H[2,2]=-integrate(dens.der,s=S,n=n,d=up.es,lower=0,upper=Inf)$value
		delta=solve(H)%*%c(rhs1,rhs2)/sqrt(n)
		low.es=low.es-delta[1];up.es=up.es-delta[2]
		if(abs(delta[1])+abs(delta[2])<eps) break
	}
	return(c(1/up.es,1/low.es))
}

set.seed(ss)
tim.st=Sys.time()
alpha=1-lambda
d0=1/kap0
kap=seq(from=.1,to=1.5,length=N) #the alternative kappas, the null is kap0, the default kap0=0.5
par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
#Equal-tailed test	
	thL.ET=qt(alpha/2,df=n-1,ncp=sqrt(n)/kap0)
	thU.ET=qt(1-alpha/2,df=n-1,ncp=sqrt(n)/kap0)
	pow.ET=1-pt(thU.ET,df=n-1,ncp=sqrt(n)/kap)+pt(thL.ET,df=n-1,ncp=sqrt(n)/kap)
	plot(kap,pow.ET,type="l",lwd=3,ylim=c(0,1),xlim=c(.1,1.6),xlab="Coefficient of variation",ylab="Power, probability")
	segments(-100,alpha,1000,alpha,col=2)
	nc=sqrt(n)/kap0
	segments(kap0,-1,kap0,2)
	text(kap[N]+.02,pow.ET[N],"ET",cex=1.5,adj=0)
	text(.7,.8,paste("n =",n),cex=2,font=3)
	
#Unbiased test		
	H=matrix(ncol=2,nrow=2)
	SL=thL.ET;SU=thU.ET
	for(it in 1:maxit)
	{
		rhs1=pt(SU,df=n-1,ncp=sqrt(n)/kap0)-pt(SL,df=n-1,ncp=sqrt(n)/kap0)-lambda
		rhs2=integrate(cdf.der1.es,s=SL,n=n,d=d0,lower=0,upper=Inf)$value-integrate(cdf.der1.es,s=SU,n=n,d=d0,lower=0,upper=Inf)$value
		H[1,1]=-dt(SL,df=n-1,ncp=sqrt(n)/kap0);H[1,2]=dt(SU,df=n-1,ncp=sqrt(n)/kap0)
		H[2,1]=integrate(f.der,s=SL,n=n,d=d0,lower=0,upper=Inf)$value
		H[2,2]=-integrate(f.der,s=SU,n=n,d=d0,lower=0,upper=Inf)$value
		delta=solve(H)%*%c(rhs1,rhs2)
		if(abs(delta[1])+abs(delta[2])<eps) break
		SL=SL-delta[1];SU=SU-delta[2]		
	}
	pow.an.UNB=1-pt(SU,df=n-1,ncp=sqrt(n)/kap)+pt(SL,df=n-1,ncp=sqrt(n)/kap)
	lines(kap,pow.an.UNB,lty=2,lwd=3)
	text(kap[N]+.02,pow.an.UNB[N],"UNB",cex=1.5,adj=0)
	SL.UNB=SL;SU.UNB=SU
		
#Density level test	

	for(it in 1:maxit)
	{
		rhs1=pt(SU,df=n-1,ncp=sqrt(n)*d0)-pt(SL,df=n-1,ncp=sqrt(n)*d0)-lambda
		rhs2=integrate(dens,s=SL,n=n,d=d0,lower=0,upper=Inf)$value-integrate(dens,s=SU,n=n,d=d0,lower=0,upper=Inf)$value
		H[1,1]=-dt(SL,df=n-1,ncp=sqrt(n)*d0);H[1,2]=dt(SU,df=n-1,ncp=sqrt(n)*d0)
		H[2,1]=integrate(dens.prime,s=SL,n=n,d=d0,lower=0,upper=Inf)$value
		H[2,2]=-integrate(dens.prime,s=SU,n=n,d=d0,lower=0,upper=Inf)$value
		delta=solve(H)%*%c(rhs1,rhs2)
		if(abs(delta[1])+abs(delta[2])<eps) break
		SL=SL-delta[1];SU=SU-delta[2]
	}	
	SL.DL=SL;SU.DL=SU
	pow.an.DL=1-pt(SU,df=n-1,ncp=sqrt(n)/kap)+pt(SL,df=n-1,ncp=sqrt(n)/kap)
	
	lines(kap,pow.an.DL,lwd=3,lty=3)
	
	legend("bottomright",c("1:Equal-tail (ET)","2:Unbiased (UNB)","3:Density level (DL)"),lty=1:3,lwd=3,cex=1.5,bg="gray96")
	text(kap[N]+.02,pow.an.DL[N],"DL",cex=1.5,adj=0)
	
#Check with simulations	based on rejection rule
	kapSIM_RR=.75;kapSIM_CI=1	
	covET=covUNB=rep(0,nSim)
	
	Y=matrix(rnorm(nSim*n,mean=1/kapSIM_RR),nrow=nSim,ncol=n)
	mY=rowSums(Y)
	mY2=rowSums(Y^2)
	sdY=sqrt((mY2-mY^2/n)/(n-1))
	S=sqrt(n)*mY/n/sdY
	powET=1-mean(S>thL.ET & S<thU.ET)
	powUNB=1-mean(S>SL.UNB & S<SU.UNB)
	powDL=1-mean(S>SL.DL & S<SU.DL)
	text(kapSIM_RR,powET,"1",cex=1.5)
	text(kapSIM_RR,powUNB,"2",cex=1.5)
	text(kapSIM_RR,powDL,"3",cex=1.5)	
	#print(c(kapSIM_RR,powET,pow_UNB,powDL))
	Y=matrix(rnorm(nSim*n,mean=1/kapSIM_CI),nrow=nSim,ncol=n)
	mY=rowSums(Y)
	mY2=rowSums(Y^2)
	sdY=sqrt((mY2-mY^2/n)/(n-1))
	S=sqrt(n)*mY/n/sdY
	for(i in 1:nSim)
	{
		ci=CI.ET(sampleES=S[i]/sqrt(n),n=n,lambda=1-alpha,maxit=100,eps=10^-7)
		covET[i]=as.numeric(ci[1]<kap0 & kap0<ci[2])
		ci=CI.unb(sampleES=S[i]/sqrt(n),n=n,lambda=1-alpha,low.unb=ci[1],up.unb=ci[2],maxit=100,eps=10^-7)
		covUNB[i]=as.numeric(ci[1]<kap0 & kap0<ci[2])		
	}
	cET=mean(covET)
	cUNB=mean(covUNB)		
	text(kapSIM_CI,1-cET,"1",cex=1.5)
	text(kapSIM_CI,1-cUNB,"2",cex=1.5)
	#print(c(kapSIM_CI,1-cET,1-cUNB))	
Sys.time()-t0
}
