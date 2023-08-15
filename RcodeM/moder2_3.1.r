moder2_3.1 <-
function(n=10,p=3,N=100)
{
dump("moder2_3.1","c:\\Projects\\Mode\\moder2_3.1.r")
#install.packages("hypergeo")

library(hypergeo)
t0=Sys.time()

dflog1=function(x,ro2,n,p)
{
	H0=Re(hypergeo((n-1)/2,(n-1)/2,p/2,x*ro2))
	H1=Re(hypergeo((n-1)/2+1,(n-1)/2+1,p/2+1,x*ro2))
	return(-(n-1)/2/(1-ro2)+x/2*(n-1)^2/p*H1/H0)		
}

fr2=function(x,ro2,n,p)
{
	fr2=gamma((n-1)/2)/gamma((n-p-1)/2)/gamma(p/2)*(1-ro2)^((n-1)/2)*(x)^((p-2)/2)*(1-x)^((n-p-3)/2)
	f0=Re(hypergeo((n-1)/2,(n-1)/2,p/2,x*ro2))
	fr2=fr2*f0 #density of r^2	
	return(fr2)
}
	
F1=function(r2,ro2,n,p)
{
	int.F1=function(x,ro2,n,p) dflog1(x,ro2,n,p)*fr2(x,ro2,n,p)
	integrate(int.F1,ro2=ro2,n=n,p=p,lower=0,upper=r2)$value
}

MCL.min=function(r2,ro2,n,p) ro2*F1(r2,ro2,n,p)

par(mfrow=c(1,2),mar=c(4,4.25,3,1),cex.main=1.5,cex.lab=1.5)
r2s=seq(from=.001,to=.9,length=200);nr2s=length(r2s)
r2.MCLa=rep(NA,nr2s)
p=c(3,7)
k=0
for(n in c(10,50))
{
	k=k+1
	for(i in 1:nr2s)
		r2.MCLa[i]=optimize(MCL.min,r2=r2s[i],n=n,p=p[k],lower=0,upper=1)$minimum
	r2.OP=1-(n-3)/(n-p[k]-1)*(1-r2s)*Re(hypergeo(1,1,(n-p[k]+1),1-r2s))	
	matplot(r2s,cbind(r2s,r2.MCLa,r2.OP),ylim=c(-.5,1),col=1,lwd=c(1,3,3),type="l",xlab="Observed MCC",ylab="Estimate of the true MCC")
	mtext(side=3,paste("n = ",n,", p = ",p[k],sep=""),cex=2,font=2,line=1)
	segments(-1,0,2,0)
	legend("topleft",c("Observed MCC","MCL estimator","Olkin-Pratt estimator"),lwd=c(1,3,3),lty=c(1,2,3),cex=1.4,bg="gray96")		
}
print(Sys.time())
}
