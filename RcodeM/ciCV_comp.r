ciCV_comp <-
function(job=2,kap0=.5,n=10,lambda=0.95,maxit=200,eps=0.00001,nSim=1000,ss=6)
{
dump(c("ciCV_comp","ciCV_comp_out"),"c:\\Projects\\Mode\\ciCV_comp.r")
#library(cvcqv)
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
		der1=-sqrt(n)/sqrt(2*pi)
		low.es=low.es-delta1
		#upper limit
		cdf0=integrate(cdf.es,s=S,n=n,d=up.es,lower=0,upper=Inf)$value
		der1=integrate(cdf.der1.es,s=S,n=n,d=up.es,lower=0,upper=Inf)$value
		der1=-sqrt(n)/sqrt(2*pi)
		delta2=(cdf0-alpha/2)/der1
		up.es=up.es-delta2
		if((abs(delta1)+abs(delta2))<eps) break
		#print(c(it,low.es,up.es))
	}
return(c(1/up.es,1/low.es))
}

CI.short=function(sampleES,n,lambda,low.SHORT,up.SHORT,maxit=100,eps=10^-7)
{
	S=sqrt(n)*sampleES
	II=0.001*matrix(c(1,0,0,1),ncol=2)
	H=matrix(ncol=2,nrow=2)
	for(it in 1:maxit)
	{
		TL=integrate(cdf.es,s=S,n=n,d=low.SHORT,lower=0,upper=Inf)$value
		TU=integrate(cdf.es,s=S,n=n,d=up.SHORT,lower=0,upper=Inf)$value
		T1L=integrate(cdf.der1.es,s=S,n=n,d=low.SHORT,lower=0,upper=Inf)$value
		T1U=integrate(cdf.der1.es,s=S,n=n,d=up.SHORT,lower=0,upper=Inf)$value
		T2L=integrate(cdf.der2.es,s=S,n=n,d=low.SHORT,lower=0,upper=Inf)$value
		T2U=integrate(cdf.der2.es,s=S,n=n,d=up.SHORT,lower=0,upper=Inf)$value
		
		rhs1=TL-TU-lambda
		rhs2=low.SHORT^2*T1L-up.SHORT^2*T1U
		#print(c(it,low.SHORT,up.SHORT,rhs1,rhs2))
		H[1,1]=T1L;H[1,2]=-T1U
		H[2,1]=2*low.SHORT*T1L+low.SHORT^2*T2L
		H[2,2]=-2*up.SHORT*T1U-up.SHORT^2*T2U
		delta=solve(H+II)%*%c(rhs1,rhs2)/sqrt(n)
		low.SHORT=low.SHORT-delta[1];up.SHORT=up.SHORT-delta[2]
		
		if(abs(rhs1)+abs(rhs2)<eps) break		
	}		
	return(c(1/up.SHORT,1/low.SHORT))
}

CI.unb=function(sampleES,n,lambda,low.unb,up.unb,maxit=100,eps=10^-7)
#the same as for effect size but returns reciprocals
#low.unb and up.unb must be low and up ET for CV (kappa)
{
	S=sqrt(n)*sampleES
	II=0.001*matrix(c(1,0,0,1),ncol=2)
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
		delta=solve(H+II)%*%c(rhs1,rhs2)/sqrt(n)
		low.es=low.es-delta[1];up.es=up.es-delta[2]
		if(abs(delta[1])+abs(delta[2])<eps) break
	}
	return(c(1/up.es,1/low.es))
}
if(job==1)
{
	set.seed(ss)
	t0=Sys.time()
	alpha=1-lambda
	d=1/kap0
	nss=seq(from=5,to=15,by=3);LA=length(nss)
	kelley=mckay=vangel=miller=ci_ET=ci_UNB=ci_SHORT=matrix(nrow=nSim,ncol=2)
	covpr=matrix(nrow=LA,ncol=6)
	for(iss in 1:LA)
	{
		n=nss[iss]
		for(isim in 1:nSim)
		{
			Y=rnorm(n,mean=d,sd=1)
			kappa.hat=sd(Y)/mean(Y)
		#	"kelley_ci", "mckay_ci", "miller_ci", "vangel_ci", "mahmoudvand_hassani_ci", "equal_tailed_ci", "shortest_length_ci", "normal_approximation_ci","norm_ci","basic_ci"
			oc=CoefVarCI$new(x=Y,alpha=alpha,R=1000,digits=4,correction=TRUE)
			
			ci=oc$mckay_ci()$statistics/100
			ocR=c(ci$lower,ci$upper)		
			mckay[isim,]=ocR			
				
			ci=oc$miller_ci()$statistics/100
			ocR=c(ci$lower,ci$upper)		
			miller[isim,]=ocR			
			
			ci=oc$vangel_ci()$statistics/100
			ocR=c(ci$lower,ci$upper)		
			vangel[isim,]=ocR	
		
			ci_ET[isim,]=CI.ET(sampleES=1/kappa.hat,n=n,lambda=1-alpha)			
			ci_UNB[isim,]=CI.unb(sampleES=1/kappa.hat,n=n,lambda=1-alpha,low.unb=ci_ET[isim,1],up.unb=ci_ET[isim,2])
			ci_SHORT[isim,]=CI.short(sampleES=1/kappa.hat,n=n,lambda=1-alpha,low.SHORT=1/ci_ET[isim,2],up.SHORT=1/ci_ET[isim,1])
		}	
			
		a5=cbind(mckay,miller,vangel,ci_ET,ci_UNB,ci_SHORT)
		a5[a5==NaN]=NA
		for(i in 1:3) covpr[iss,i]=mean(a5[,1+(i-1)*2]<kap0 & a5[,2+(i-1)*2]>kap0,na.rm=T)
		for(i in 4:6) 
		{
			kL=a5[,1+(i-1)*2];kU=a5[,2+(i-1)*2]
			covpr[iss,i]=mean(kL<kap0 & kap0<kU)+mean(kap0>kL & kL>kU)+mean(kap0<kU & kL>kU)
		}	
	}
	print(Sys.time()-t0)
	out=as.data.frame(covpr,row.names=nss)
	names(out)=c("McKay","Miller","Vangel","ET","UNB","SHORT")
	return(out)
	#ciCV_comp_out=ciCV_comp(nSim=1000)	
}
if(job==2)
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	n=as.numeric(row.names(ciCV_comp_out))
	covpr=as.matrix(ciCV_comp_out)
	matplot(n,covpr,type="o",col=1,lwd=2,ylim=c(0.75,1),lty=1,xlab="Sample size, n",ylab="Coverage probablity")
	segments(0,lambda,20,lambda,col=2)
	legend("bottomright",paste(as.character(1:6),names(ciCV_comp_out)),cex=1.5,bg="gray96")
	text(10,.98,paste("True CV=",kap0,", confidence level=",lambda*100,"%",sep=""),cex=1.5,font=3)
}
}
ciCV_comp_out <-
structure(list(McKay = c(0.77629826897470044, 0.90723822909346452, 
0.93544857768052514, 0.94799999999999995), Miller = c(0.82650000000000001, 
0.86950000000000005, 0.88700000000000001, 0.91300000000000003
), Vangel = c(0.82346491228070173, 0.92332065906210392, 0.94207477619799895, 
0.94850000000000001), ET = c(0.94650000000000001, 0.95450000000000002, 
0.95499999999999996, 0.95450000000000002), UNB = c(0.94600000000000006, 
0.94850000000000001, 0.94750000000000001, 0.95399999999999996
), SHORT = c(0.95100000000000007, 0.95649999999999991, 0.95599999999999996, 
0.96399999999999997)), class = "data.frame", row.names = c("5", 
"8", "11", "14"))
