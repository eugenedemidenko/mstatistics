testES <-
function(job=1,n=10,ES.null=0.5,ES.alt=1,sampleES=0.25,alpha=0.05,maxit=10,eps=10^-6)
{
dump("testES","c:\\Projects\\Mode\\testES.r")

cdf.der1.es=function(w,s,n,d)	-2*sqrt(n)*w*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
f.der=function(w,s,n,d)	2*sqrt(n)*w^2/sqrt(n-1)*(s*w/sqrt(n-1)-d*sqrt(n))*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
dens=function(w,s,n,d)	2/sqrt(n-1)*w^2*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
dens.prime=function(w,s,n,d) -2/(n-1)*w^3*(s*w/sqrt(n-1)-sqrt(n)*d)*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)



if(job==1)#powers
{

#Equal-tailed test	
	SL=SL.ET=qt(alpha/2,df=n-1,ncp=sqrt(n)*ES.null)
	SU=SU.ET=qt(1-alpha/2,df=n-1,ncp=sqrt(n)*ES.null)
	pow.ET=1-pt(SU,df=n-1,ncp=sqrt(n)*ES.alt)+pt(SL,df=n-1,ncp=sqrt(n)*ES.alt)
	
#Unbiased test		
	H=matrix(ncol=2,nrow=2)
	for(it in 1:maxit)
	{
		rhs1=pt(SU,df=n-1,ncp=sqrt(n)*ES.null)-pt(SL,df=n-1,ncp=sqrt(n)*ES.null)-(1-alpha)
		rhs2=integrate(cdf.der1.es,s=SL,n=n,d=ES.null,lower=0,upper=Inf)$value-integrate(cdf.der1.es,s=SU,n=n,d=ES.null,lower=0,upper=Inf)$value
		H[1,1]=-dt(SL,df=n-1,ncp=sqrt(n)*ES.null);H[1,2]=dt(SU,df=n-1,ncp=sqrt(n)*ES.null)
		H[2,1]=integrate(f.der,s=SL,n=n,d=ES.null,lower=0,upper=Inf)$value
		H[2,2]=-integrate(f.der,s=SU,n=n,d=ES.null,lower=0,upper=Inf)$value
		delta=solve(H)%*%c(rhs1,rhs2)
		if(abs(delta[1])+abs(delta[2])<eps) break
		SL=SL-delta[1];SU=SU-delta[2]		
	}
	pow.UNB=1-pt(SU,df=n-1,ncp=sqrt(n)*ES.alt)+pt(SL,df=n-1,ncp=sqrt(n)*ES.alt)			
#Density level test	
	for(it in 1:maxit)
	{
		rhs1=pt(SU,df=n-1,ncp=sqrt(n)*ES.null)-pt(SL,df=n-1,ncp=sqrt(n)*ES.null)-(1-alpha)
		rhs2=integrate(dens,s=SL,n=n,d=ES.null,lower=0,upper=Inf)$value-integrate(dens,s=SU,n=n,d=ES.null,lower=0,upper=Inf)$value
		H[1,1]=-dt(SL,df=n-1,ncp=sqrt(n)*ES.null);H[1,2]=dt(SU,df=n-1,ncp=sqrt(n)*ES.null)
		H[2,1]=integrate(dens.prime,s=SL,n=n,d=ES.null,lower=0,upper=Inf)$value
		H[2,2]=-integrate(dens.prime,s=SU,n=n,d=ES.null,lower=0,upper=Inf)$value
		delta=solve(H)%*%c(rhs1,rhs2)
		if(abs(delta[1])+abs(delta[2])<eps) break
		SL=SL-delta[1];SU=SU-delta[2]
	}	
	pow.DL=1-pt(SU,df=n-1,ncp=sqrt(n)*ES.alt)+pt(SL,df=n-1,ncp=sqrt(n)*ES.alt)
	cat("ES.null=",ES.null,", ES.alt=",ES.alt," n=",n,sep="","\n")	
	cat("Power equal-tail test =",round(pow.ET,4),"\n")
	cat("Power inbiased test =",round(pow.UNB,4),"\n")
	cat("Power density level test =",round(pow.DL,4),"\n")

}
if(job==2)#p-values
{
	cat("sample ES=",sampleES,", ES.null=",ES.null," n=",n,sep="","\n")	
	S=sqrt(n)*sampleES
#Equal-tail	
	pL=pt(S,df=n-1,ncp=sqrt(n)*ES.null)
	pU=pt(S,df=n-1,ncp=sqrt(n)*ES.null,lower.tail=F)
	pET=2*min(pL,pU)
#Unbiased test			
	SL=SL.ET=qt(alpha/2,df=n-1,ncp=sqrt(n)*ES.null)
	SU=SU.ET=qt(1-alpha/2,df=n-1,ncp=sqrt(n)*ES.null)
	H=matrix(ncol=2,nrow=2)
	for(it in 1:maxit)
	{
		rhs1=pt(SU,df=n-1,ncp=sqrt(n)*ES.null)-pt(SL,df=n-1,ncp=sqrt(n)*ES.null)-(1-alpha)
		rhs2=integrate(cdf.der1.es,s=SL,n=n,d=ES.null,lower=0,upper=Inf)$value-integrate(cdf.der1.es,s=SU,n=n,d=ES.null,lower=0,upper=Inf)$value
		H[1,1]=-dt(SL,df=n-1,ncp=sqrt(n)*ES.null);H[1,2]=dt(SU,df=n-1,ncp=sqrt(n)*ES.null)
		H[2,1]=integrate(f.der,s=SL,n=n,d=ES.null,lower=0,upper=Inf)$value
		H[2,2]=-integrate(f.der,s=SU,n=n,d=ES.null,lower=0,upper=Inf)$value
		delta=solve(H)%*%c(rhs1,rhs2)
		if(abs(delta[1])+abs(delta[2])<eps) break
		SL=SL-delta[1];SU=SU-delta[2]		
	}
	FS0=pt(S,df=n-1,ncp=sqrt(n)*ES.null)
	FqL=pt(SL,df=n-1,ncp=sqrt(n)*ES.null)
	if(FS0<FqL/alpha) pUNB=alpha/FqL*FS0 else pUNB=alpha/(alpha-FqL)*(1-FS0)
#Density level test	
	for(it in 1:maxit)
	{
		rhs1=pt(SU,df=n-1,ncp=sqrt(n)*ES.null)-pt(SL,df=n-1,ncp=sqrt(n)*ES.null)-(1-alpha)
		rhs2=integrate(dens,s=SL,n=n,d=ES.null,lower=0,upper=Inf)$value-integrate(dens,s=SU,n=n,d=ES.null,lower=0,upper=Inf)$value
		H[1,1]=-dt(SL,df=n-1,ncp=sqrt(n)*ES.null);H[1,2]=dt(SU,df=n-1,ncp=sqrt(n)*ES.null)
		H[2,1]=integrate(dens.prime,s=SL,n=n,d=ES.null,lower=0,upper=Inf)$value
		H[2,2]=-integrate(dens.prime,s=SU,n=n,d=ES.null,lower=0,upper=Inf)$value
		delta=solve(H)%*%c(rhs1,rhs2)
		if(abs(delta[1])+abs(delta[2])<eps) break
		SL=SL-delta[1];SU=SU-delta[2]
	}	
	FqL=pt(SL,df=n-1,ncp=sqrt(n)*ES.null)
	if(FS0<FqL/alpha) pDL=alpha/FqL*FS0 else pDL=alpha/(alpha-FqL)*(1-FS0)
	cat("P-value equal-tail test =",round(pET,4),"\n")
	cat("Power inbiased test =",round(pUNB,4),"\n")
	cat("Power density level test =",round(pDL,4),"\n")
	

}
}
