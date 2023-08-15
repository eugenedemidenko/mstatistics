test_r <-
function(n=10,ro0=.25,ro.alt=.5,r=.1,alpha=0.05,maxit=10,ss=3,nSim=1000000)
{
dump("test_r","c:\\Projects\\Mode\\test_r.r")
#install.packages("hypergeo")
library(hypergeo)
set.seed(ss)
t0=Sys.time()

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


out.typeI=as.data.frame(matrix(ncol=5,nrow=3))
row.names(out.typeI)=c("Exact equal-tail","Fisher approx","Unbiased exact")
names(out.typeI)=c("Sim type I","Theor type I","Sim power","Theor power","p-value power")

#Simulation with ro=ro0	
X=matrix(rnorm(nSim*n),ncol=n);xm=rowMeans(X)
X2=rowSums(X^2)
U=matrix(rnorm(nSim*n),ncol=n)
Y=ro0*X+U*sqrt(1-ro0^2)
Y2=rowSums(Y^2);XY=rowSums(X*Y)
ym=rowMeans(Y)
r.obs=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))

#Simulation with ro=ro.alt	
X=matrix(rnorm(nSim*n),ncol=n);xm=rowMeans(X)
X2=rowSums(X^2)
U=matrix(rnorm(nSim*n),ncol=n)
Y=ro.alt*X+U*sqrt(1-ro.alt^2)
Y2=rowSums(Y^2);XY=rowSums(X*Y)
ym=rowMeans(Y)
r.alt=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))

#Fisher approx
qL.Z=Z.inv(Z(ro0)+qnorm(alpha/2)/sqrt(n-3))
qU.Z=Z.inv(Z(ro0)+qnorm(1-alpha/2)/sqrt(n-3))
out.typeI[2,1]=1-mean(r.obs>qL.Z & r.obs<qU.Z)	
out.typeI[2,2]=1-(cdf.ro(r=qU.Z,ro=ro0,n=n)-cdf.ro(r=qL.Z,ro=ro0,n=n))
out.typeI[2,3]=1-mean(r.alt>qL.Z & r.alt<qU.Z)	
out.typeI[2,4]=1+cdf.ro(r=qL.Z,ro=ro.alt,n=n)-cdf.ro(r=qU.Z,ro=ro.alt,n=n)
pv.Z=2*pnorm(-abs(Z(r.alt)-Z(ro0))*sqrt(n-3))
out.typeI[2,5]=mean(pv.Z<alpha)

#Exact equal-tail	
qL.ET=qL.Z;qU.ET=qU.Z
for(iss in 0:maxit)
{
	qL.ET=qL.ET-(cdf.ro(r=qL.ET,ro=ro0,n=n)-alpha/2)/dens.ro(r=qL.ET,ro=ro0,n=n)
	qU.ET=qU.ET-(cdf.ro(r=qU.ET,ro=ro0,n=n)-(1-alpha/2))/dens.ro(r=qU.ET,ro=ro0,n=n)		
}	
out.typeI[1,1]=1-mean(r.obs>qL.ET & r.obs<qU.ET)
out.typeI[1,2]=1-(cdf.ro(r=qU.ET,ro=ro0,n=n)-cdf.ro(r=qL.ET,ro=ro0,n=n))
out.typeI[1,3]=1-mean(r.alt>qL.ET & r.alt<qU.ET)	
out.typeI[1,4]=1+cdf.ro(r=qL.ET,ro=ro.alt,n=n)-cdf.ro(r=qU.ET,ro=ro.alt,n=n)
Fr=rep(NA,nSim)
for(i in 1:nSim)
	Fr[i]=cdf.ro(r=r.alt[i],ro=ro0,n=n)

pv.ET=2*pmin(Fr,1-Fr)
out.typeI[1,5]=mean(pv.ET<alpha)

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
out.typeI[3,1]=1-mean(r.obs>qL.UNB & r.obs<qU.UNB)
out.typeI[3,2]=1-(cdf.ro(r=qU.UNB,ro=ro0,n=n)-cdf.ro(r=qL.UNB,ro=ro0,n=n))
out.typeI[3,3]=1-mean(r.alt>qL.UNB & r.alt<qU.UNB)	
out.typeI[3,4]=1+cdf.ro(r=qL.UNB,ro=ro.alt,n=n)-cdf.ro(r=qU.UNB,ro=ro.alt,n=n)
Frq=cdf.ro(r=qL.UNB,ro=ro0,n=n)
ind=(Fr <= Frq/alpha)
pv.UNB=ind*alpha/Frq*Fr+(1-ind)*alpha/(alpha-Frq)*(1-Fr)
out.typeI[3,5]=mean(pv.UNB<alpha)
print(Sys.time()-t0)	
out.typeI

}
