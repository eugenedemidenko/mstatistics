test_rpow <-
function(n=10,ro0=.8,ro.alt=.9,alpha=0.05,maxit=10)
{
dump("test_rpow","c:\\Projects\\Mode\\test_rpow.r")
#install.packages("hypergeo")
library(hypergeo)

Z=function(x) 0.5*log((1+x)/(1-x))
Z.inv=function(z) (exp(2*z)-1)/(exp(2*z)+1)

dens.ro=function(r,ro,n) #r density
{
	COF.num=(n-2)*gamma(n-1)*(1-ro^2)^(n/2-.5)*(1-r^2)^(n/2-2)
	COF.den=sqrt(2*pi)*gamma(n-0.5)*(1-ro*r)^(n-1.5)
	fr=COF.num/COF.den*hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5)
	return(Re(fr))
}
cdf.ro=function(r,ro,n) integrate(f=dens.ro,ro=ro,n=n,lower=-1,upper=r)$value

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
cdf.dero=function(r,ro,n) integrate(f=densder.ro,ro=ro,n=n,lower=-1,upper=r)$value

par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
ros=seq(from=0.5,to=0.95,by=0.01);nros=length(ros)
 
#Exact equal-tail	
qL.ET=qU.ET=ro0
for(iss in 0:maxit)
{
	qL.ET=qL.ET-(cdf.ro(r=qL.ET,ro=ro0,n=n)-alpha/2)/dens.ro(r=qL.ET,ro=ro0,n=n)
	qU.ET=qU.ET-(cdf.ro(r=qU.ET,ro=ro0,n=n)-(1-alpha/2))/dens.ro(r=qU.ET,ro=ro0,n=n)	
	#print(c(iss,qL,qU))	
}	

#Fisher approx
qL.Z=Z.inv(Z(ro0)+qnorm(alpha/2)/sqrt(n-3))
qU.Z=Z.inv(Z(ro0)+qnorm(1-alpha/2)/sqrt(n-3))


#Unbiased exact test
#Unbiased test q1=qL and q2=qU starting from Fisher-Z
q1=qL.ET;q2=qU.ET
H=matrix(ncol=2,nrow=2)
for(it in 1:maxit)
{
	rh1=cdf.ro(q2,ro=ro0,n=n)-cdf.ro(q1,ro=ro0,n=n)-(1-alpha)
	rh2=cdf.dero(q1,ro=ro0,n=n)-cdf.dero(q2,ro=ro0,n=n)
	H[1,1]=-dens.ro(q1,ro=ro0,n=n);H[1,2]=dens.ro(q2,ro=ro0,n=n)
	H[2,1]=densder.ro(q1,ro=ro0,n=n);H[2,2]=-densder.ro(q2,ro=ro0,n=n)
	delta=solve(H)%*%c(rh1,rh2)
	q12=c(q1,q2)-delta
	if(max(abs(delta))<10^-7) break
	q1=q12[1];q2=q12[2]	
}
qL.UNB=q1;qU.UNB=q2

powET=powZ=powUNB=rep(NA,nros)
for(i in 1:nros)
{
	powET[i]=1+cdf.ro(r=qL.ET,ro=ros[i],n=n)-cdf.ro(r=qU.ET,ro=ros[i],n=n)
	powZ[i]=1+cdf.ro(r=qL.Z,ro=ros[i],n=n)-cdf.ro(r=qU.Z,ro=ros[i],n=n)
	powUNB[i]=1+cdf.ro(r=qL.UNB,ro=ros[i],n=n)-cdf.ro(r=qU.UNB,ro=ros[i],n=n)
}
matplot(ros,cbind(powET,powZ,powUNB),ylim=c(0,.4),col=1,lwd=2,type="l",xlab="Alternative ro",ylab="Power")
segments(-1,alpha,1,alpha,lty=2)
segments(ro0,-1,ro0,1,lty=2)
text(rep(.496,3),c(powET[1],powZ[1],powUNB[1])-.005,as.character(1:3),pos=3,cex=1.4,font=2)
legend("bottomleft",c("Equal-tail","Fisher Z","Unbiased"),cex=1.5,lty=1:3,lwd=2,bg="grey96")
#cbind(ros,powET,powZ,powUNB)

}
