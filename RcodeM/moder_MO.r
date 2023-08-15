moder_MO <-
function(n=10)
{
dump("moder_MO","c:\\Projects\\Mode\\moder_MO.r")
#install.packages("hypergeo")

library(hypergeo)

lndens.ro=function(ro,r,n) #log r density up to a const
{
	fr=(n-1)/2*log(1-ro^2)-(2*n-3)/2*log(1-ro*r)+log(hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5))
	return(Re(fr))
}
cdf.dero=function(r,ro,n)  #derivative F'
{
	cdf.val=integrate(f=densder.ro,ro=ro,n=n,lower=-1,upper=r)$value
	return(cdf.val)	
}
densder.ro=function(r,ro,n)
{
	COF.num=(n-2)*gamma(n-1)*(1-ro^2)^(n/2-.5)*(1-r^2)^(n/2-2)
	COF.den=sqrt(2*pi)*gamma(n-0.5)*(1-ro*r)^(n-1.5)
	f0=hypergeo(0.5,0.5,n-0.5,r*ro/2+0.5)
	fr=Re(COF.num/COF.den*f0)
			
	f1=hypergeo(0.5+1,0.5+1,n-0.5+1,r*ro/2+0.5)
	hr=-ro*(n-1)/(1-ro^2)+r*(2*n-3)/2/(1-ro*r)+r/4/(2*n-1)*Re(f1/f0)
	return(fr*hr)		
}

#MC, MO and Olkin-Pratt estimator

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
