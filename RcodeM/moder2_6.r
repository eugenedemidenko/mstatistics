moder2_6 <-
function(N=100,lambda=0.95,PEN=1000)
{
dump("moder2_6","c:\\Projects\\Mode\\moder2_6.r")

#install.packages("hypergeo")
#install.packages("MBESS")
library(MBESS)
library(hypergeo)
t0=Sys.time()

fr2=function(x,ro2,n,p)
{
	fr2=gamma((n-1)/2)/gamma((n-p-1)/2)/gamma(p/2)*(1-ro2)^((n-1)/2)*(x)^((p-2)/2)*(1-x)^((n-p-3)/2)
	f0=Re(hypergeo((n-1)/2,(n-1)/2,p/2,x*ro2))
	fr2=fr2*f0 #pdf of r2	
	return(fr2)
}
	
F1=function(r2,ro2,n,p)  #cdf of r2	
{
	int.F1=function(x,ro2,n,p) dflog1(x,ro2,n,p)*fr2(x,ro2,n,p)
	integrate(int.F1,ro2=ro2,n=n,p=p,lower=0,upper=r2)$value
}
dflog1=function(x,ro2,n,p) #1st der of pdf
{
	H0=Re(hypergeo((n-1)/2,(n-1)/2,p/2,x*ro2))
	H1=Re(hypergeo((n-1)/2+1,(n-1)/2+1,p/2+1,x*ro2))
	return(-(n-1)/2/(1-ro2)+x/2*(n-1)^2/p*H1/H0)		
}

MCL.min=function(r2,ro2,n,p) ro2*F1(r2,ro2,n,p) #function to minimize to get MCL estimate of r2

fnl=function(q,r2,n,p,lambda=0.95,PEN=1000) #penalyzed function to minimize for the CI
{
	eq1=integrate(fr2,ro2=q[1],n,p,lower=0.0001,upper=r2)$value-integrate(fr2,ro2=q[2],n,p,lower=0.0001,upper=r2)$value-lambda
	return(c(PEN*eq1^2+log(q[2]/q[1])^2))	
}

par(mfrow=c(1,2),cex.lab=1.4,cex.main=1.4)
p=c(2,4);n=c(10,30)
r2s=r2.MCL=seq(from=0.05,to=0.95,length=N)
x=seq(from=0,to=1,length=N)
qLU.MBESS=qLU.PEN=matrix(ncol=2,nrow=N)
for(ic in 1:2)
{
	for(i in 1:N)
	{
		out=ci.R2(R2=r2s[i],df.1=p[ic],df.2=n[ic]-p[ic]-1,Random.Regressors=TRUE,conf.level=lambda)
		qLU.MBESS[i,]=c(qL=out$Lower.Conf.Limit.R,qU=out$Upper.Conf.Limit.R) #MBESS
	
		ci=optim(par=c(0.1,0.9),fnl,lower=0.0001,upper=1,method="L-BFGS-B",r2=r2s[i],n=n[ic],p=p[ic],lambda=lambda,PEN=PEN)$par	#penalyzed function for CI on the log scale 
		qLU.PEN[i,1]=min(ci);qLU.PEN[i,2]=max(ci)		
		r2.MCL[i]=optimize(MCL.min,r2=r2s[i],n=n[ic],p=p[ic],lower=0,upper=1)$minimum #MCL estimate
			
	}
	matplot(r2s,cbind(qLU.PEN,qLU.MBESS),col=1,lty=c(1,1,2,2),ylim=c(0,1),xlim=c(0,1),lwd=2,type="l",xlab="Observed MSCC",ylab="Confidence interval")
	lines(r2s,r2.MCL,col=2)
	title(paste("n = ",n[ic],", p = ",p[ic],sep=""))
	segments(-1,-1,2,2,lty=4)
	legend("topleft",c("Short exact CI","Approximate CI","MCL estimate"),lty=c(1:2,1),col=c(1,1,2),lwd=c(2,2,1),cex=1.5,bg="gray93")		
}
print(Sys.time()-t0)
}
