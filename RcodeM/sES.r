sES <-
function(n=10,sampleES=1.5,alpha=0.05,maxit=100,eps=10^-6)
{
dump("sES","c:\\Projects\\Mode\\sES.r")
#install.packages("MBESS")
#If you get the message 'package MBESS is not abvailable' follow the these steps:
#install.packages("remotes")
#library(remotes)
#install_github("cran/MBESS")
#If you get an error messages again install a newer version of R
library(MBESS)

cdf.es=function(w,s,n,d)2*w*pnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
cdf.der1.es=function(w,s,n,d)-2*sqrt(n)*w*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
cdf.der2.es=function(w,s,n,d)-2*n*w*(s*w/sqrt(n-1)-d*sqrt(n))*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
cdf.der3.es=function(w,s,n,d)	-2*n^1.5*w*((s*w/sqrt(n-1)-d*sqrt(n))^2-1)*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
dens=function(w,s,n,d)2/sqrt(n-1)*w^2*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
dens.der=function(w,s,n,d)2*sqrt(n)/sqrt(n-1)*w^2*(s*w/sqrt(n-1)-sqrt(n)*d)*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
dens.prime=function(w,s,n,d) -2/(n-1)*w^3*(s*w/sqrt(n-1)-sqrt(n)*d)*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)
dens.der.prime=function(w,s,n,d) 2*sqrt(n)/(n-1)*w^3*(1+(s*w/sqrt(n-1)-sqrt(n)*d)^2)*dnorm(s*w/sqrt(n-1)-d*sqrt(n))*dchisq(w^2,df=n-1)

#------------------------------------------------------MBESS CI
outMBESS=ci.sm(sm=sampleES,N=n,conf.level=1-alpha,message=F)
MBESS_low=outMBESS$Lower.Conf.Limit.Standardized.Mean
MBESS_up=outMBESS$Upper.Conf.Limit.Standardized.Mean		

S=sqrt(n)*sampleES

#-----------------------------------------------------Equal-tail
# Johnson and Welch (1940) approximations
low.es=(S-qnorm(1-alpha/2)*sqrt(1+S^2/2/n))/sqrt(n)
up.es=(S-qnorm(alpha/2)*sqrt(1+S^2/2/n))/sqrt(n)

#print("ET")
for(it in 1:maxit)
{
	#print(c(it,low.es,up.es))
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
}

#------------------------------------------------------short CI
#print("SHORT")
H=matrix(ncol=2,nrow=2)
low.sh=low.es;up.sh=up.es
for(it in 1:maxit)
{
	#print(c(it,low.sh,up.sh))
	rhs1=integrate(cdf.es,s=S,n=n,d=low.sh,lower=0,upper=Inf)$value-integrate(cdf.es,s=S,n=n,d=up.sh,lower=0,upper=Inf)$value - (1-alpha)
	H[1,1]=integrate(cdf.der1.es,s=S,n=n,d=low.sh,lower=0,upper=Inf)$value
	H[1,2]=-integrate(cdf.der1.es,s=S,n=n,d=up.sh,lower=0,upper=Inf)$value
	rhs2=H[1,1]+H[1,2]
	H[2,1]=integrate(cdf.der2.es,s=S,n=n,d=low.sh,lower=0,upper=Inf)$value
	H[2,2]=-integrate(cdf.der2.es,s=S,n=n,d=up.sh,lower=0,upper=Inf)$value
	delta=solve(H)%*%c(rhs1,rhs2)/sqrt(n)
	low.sh=low.sh-delta[1];up.sh=up.sh-delta[2]
	if(abs(delta[1])+abs(delta[2])<eps) break
}

#-----------------------------------------------------unbiased CI
#print("UNB")
H=matrix(ncol=2,nrow=2)
low.unb=low.es;up.unb=up.es
for(it in 1:maxit)
{
	#print(c(it,low.unb,up.unb))
	rhs1=pt(S,df=n-1,ncp=sqrt(n)*low.unb)-pt(S,df=n-1,ncp=sqrt(n)*up.unb) - (1-alpha)
	H[1,1]=integrate(cdf.der1.es,s=S,n=n,d=low.unb,lower=0,upper=Inf)$value
	H[1,2]=-integrate(cdf.der1.es,s=S,n=n,d=up.unb,lower=0,upper=Inf)$value
	rhs2=integrate(dens,s=S,n=n,d=low.unb,lower=0,upper=Inf)$value-integrate(dens,s=S,n=n,d=up.unb,lower=0,upper=Inf)$value
	H[2,1]=integrate(dens.der,s=S,n=n,d=low.unb,lower=0,upper=Inf)$value
	H[2,2]=-integrate(dens.der,s=S,n=n,d=up.unb,lower=0,upper=Inf)$value
	delta=solve(H)%*%c(rhs1,rhs2)/sqrt(n)
	low.unb=low.unb-delta[1];up.unb=up.unb-delta[2]
	if(abs(delta[1])+abs(delta[2])<eps) break
}
out=as.data.frame(cbind(c(MBESS_low,MBESS_up),c(low.es,up.es),c(low.sh,up.sh),c(low.unb,up.unb)))
names(out)=c("MBESS","Equal-tail CI","Short CI","Unbiased CI")
row.names(out)=c("Lower","Upper")
cat("Sample ES=",sampleES,", n=",n,"\n",sep="")

#--------------------------------------------------MC estimate
MC.ES=sampleES
for(it in 1:maxit)
{
	num=integrate(cdf.der2.es,s=S,n=n,d=MC.ES,lower=0,upper=Inf)$value
	den=integrate(cdf.der3.es,s=S,n=n,d=MC.ES,lower=0,upper=Inf)$value
	delta=num/den
	MC.ES=MC.ES-delta
	#print(c(it,MC.ES,num))
	if(abs(delta)<eps) break			
}
#--------------------------------------------------MO estimate	
MO.ES=sampleES
for(it in 1:maxit)
{
	num=integrate(dens.prime,s=S,n=n,d=MO.ES,lower=0,upper=Inf)$value
	den=integrate(dens.der.prime,s=S,n=n,d=MO.ES,lower=0,upper=Inf)$value
	delta=num/den
	MO.ES=MO.ES-delta
	#print(c(it,MO.ES,num))
	if(abs(delta)<eps) break							
}	
cat("MC estimate=",MC.ES,", MO estimate=",MO.ES,"\n95% confidence limits\n",sep="")
out

}
