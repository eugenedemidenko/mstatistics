moder_MC <-
function(n=10,N=100)
{
dump("moder_MC","c:\\Projects\\Mode\\moder_MC.r")
#install.packages("hypergeo")

library(hypergeo)

dens.ro=function(r,ro,n) #r density
{
	COF.num=(n-2)*gamma(n-1)*(1-ro^2)^(n/2-.5)*(1-r^2)^(n/2-2)
	COF.den=sqrt(2*pi)*gamma(n-0.5)*(1-ro*r)^(n-1.5)
	fr=COF.num/COF.den*hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5)
	return(Re(fr))
}
der1=function(r,ro,n) #derivative of log density
{
	H0=hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5)
	H1=hypergeo(0.5+1,0.5+1,n-0.5+1,ro*r/2+0.5)
	der1=-ro*(n-1)/(1-ro^2)+r*(2*n-3)/2/(1-ro*r)+r/4/(2*n-1)*H1/H0
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
num.ro=function(r,ro,n) (der2(r,ro,n)+(der1(r,ro,n))^2)*dens.ro(r,ro,n)
num.int=function(r,ro,n) integrate(num.ro,ro=ro,n=n,lower=-1,upper=r)$value

#MC estimator illustration

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
