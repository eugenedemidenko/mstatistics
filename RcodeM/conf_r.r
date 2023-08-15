conf_r <-
function(n=10,ro0=.25,alpha=0.05,eps=10^-9,maxit=100,ss=3,nSim=1000000)
{
dump("conf_r","c:\\Projects\\Mode\\conf_r.r")
#install.packages("hypergeo")
#library(hypergeo)
install.packages("gsl")
library(gsl) #an alternative to hypergeo
hypergeo=function(a,b,cc,x)
hyperg_2F1(a,b,cc,x)

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
logdensder.ro=function(r,ro,n)
{
	f0=hypergeo(0.5,0.5,n-0.5,r*ro/2+0.5)
	f1=hypergeo(0.5+1,0.5+1,n-0.5+1,r*ro/2+0.5)
	hr=-ro*(n-1)/(1-ro^2)+r*(2*n-3)/2/(1-ro*r)+r/4/(2*n-1)*Re(f1/f0)
	return(hr)		
}
logdensder2.ro=function(r,ro,n)
{
	H0=hypergeo(0.5,0.5,n-0.5,ro*r/2+0.5)
	H1=hypergeo(0.5+1,0.5+1,n-0.5+1,ro*r/2+0.5)
	H2=hypergeo(0.5+2,0.5+2,n-0.5+2,ro*r/2+0.5)
	der2=-(1+ro^2)*(n-1)/(1-ro^2)^2+r^2*(2*n-3)/2/(1-ro*r)^2+r^2/16/(2*n-1)*(9/(2*n+1)*H2/H0-1/(2*n-1)*H1^2/H0^2)
	return(Re(der2))
}
num.ro=function(r,ro,n) (logdensder2.ro(r,ro,n)+(logdensder.ro(r,ro,n))^2)*dens.ro(r,ro,n)
num.int=function(r,ro,n) integrate(num.ro,ro=ro,n=n,lower=-1,upper=r)$value


#rs=seq(from=-.99,to=.99,by=.1);nrs=length(rs)
#print(nrs)
#par(mfrow=c(4,5),mar=c(1,1,1,1))
#ros=seq(from=-.99,to=.99,by=.01);nros=length(ros)
#F2=rep(NA,nros)
#for(j in 1:nrs)
#{
#	for(i in 1:nros) F2[i]=num.int(r=rs[j],ro=ros[i],n)
#	plot(ros,F2,col=1,lty=1,type="l")
#	segments(-1,0,1,0,col=2)
#}
#return()


out.typeI=as.data.frame(matrix(ncol=1,nrow=4))
row.names(out.typeI)=c("Fisher approx","Exact equal-tail","Unbiased exact","Short exact")
names(out.typeI)="Coverage probability" 

#Simulation with ro=ro0	
X=matrix(rnorm(nSim*n),ncol=n);xm=rowMeans(X)
X2=rowSums(X^2)
U=matrix(rnorm(nSim*n),ncol=n)
Y=ro0*X+U*sqrt(1-ro0^2)
Y2=rowSums(Y^2);XY=rowSums(X*Y)
ym=rowMeans(Y)
r.obs=(XY-n*xm*ym)/sqrt((X2-n*xm^2)*(Y2-n*ym^2))


#Fisher approx CI
roL.Z=Z.inv(Z(r.obs)-qnorm(1-alpha/2)/sqrt(n-3))
roU.Z=Z.inv(Z(r.obs)+qnorm(1-alpha/2)/sqrt(n-3))
out.typeI[1,1]=mean(roL.Z<ro0 & ro0<roU.Z)	

#Exact equal-tail CI	
roL.ET=roL.Z;roU.ET=roU.Z
for(i in 1:nSim)
{
	for(iss in 0:maxit)
	{
		deltaU=(cdf.ro(r=r.obs[i],ro=roU.ET[i],n=n)-alpha/2)/cdf.dero(r=r.obs[i],ro=roU.ET[i],n=n)
		roU.ET[i]=roU.ET[i]-deltaU
		deltaL=(cdf.ro(r=r.obs[i],ro=roL.ET[i],n=n)-(1-alpha/2))/cdf.dero(r=r.obs[i],ro=roL.ET[i],n=n)				
		roL.ET[i]=roL.ET[i]-deltaL
		if(max(deltaU,deltaL)<eps) break
	}	
}
out.typeI[2,1]=mean(roL.ET<ro0 & ro0<roU.ET)	

#Unbiased exact CI
roL.UNB=roL.ET;roU.UNB=roU.ET
H=matrix(ncol=2,nrow=2)
for(i in 1:nSim)
{
	for(iss in 0:maxit)
	{
		
		curit=iss
		eq1=cdf.ro(r=r.obs[i],ro=roL.UNB[i],n=n)-cdf.ro(r=r.obs[i],ro=roU.UNB[i],n=n)-(1-alpha)
		eq2=log(dens.ro(r=r.obs[i],ro=roL.UNB[i],n=n))-log(dens.ro(r=r.obs[i],ro=roU.UNB[i],n=n))
		H[1,1]=cdf.dero(r=r.obs[i],ro=roL.UNB[i],n=n)
		H[1,2]=-cdf.dero(r=r.obs[i],ro=roU.UNB[i],n=n)
		H[2,1]=logdensder.ro(r=r.obs[i],ro=roL.UNB[i],n=n)
		H[2,2]=-logdensder.ro(r=r.obs[i],ro=roU.UNB[i],n=n)
		iH=solve(H)
		delta=iH%*%c(eq1,eq2)
		roL.UNB[i]=roL.UNB[i]-delta[1]
		roU.UNB[i]=roU.UNB[i]-delta[2]		
		#print(c(iss,roL.UNB[i],roU.UNB[i],eq1,eq2))
		if(max(abs(delta))<eps) break
	}	
	if(curit==maxit) roL.UNB[i]=roU.UNB[i]=NA
}
out.typeI[3,1]=mean(roL.UNB<ro0 & ro0<roU.UNB,na.rm=T)	

#Short CI
roL.SHORT=roL.ET;roU.SHORT=roU.ET

for(i in 1:nSim)
{
	for(iss in 0:maxit)
	{
		curit=iss
		eq1=cdf.ro(r=r.obs[i],ro=roL.SHORT[i],n=n)-cdf.ro(r=r.obs[i],ro=roU.SHORT[i],n=n)-(1-alpha)
		eq2=cdf.dero(r=r.obs[i],ro=roL.SHORT[i],n=n)-cdf.dero(r=r.obs[i],ro=roU.SHORT[i],n=n)
		H[1,1]=cdf.dero(r=r.obs[i],ro=roL.SHORT[i],n=n)
		H[1,2]=-cdf.dero(r=r.obs[i],ro=roU.SHORT[i],n=n)
		H[2,1]=num.int(r=r.obs[i],ro=roL.SHORT[i],n=n)
		H[2,2]=-num.int(r=r.obs[i],ro=roU.SHORT[i],n=n)
		iH=solve(H)
		delta=iH%*%c(eq1,eq2)
		roL.SHORT[i]=roL.SHORT[i]-delta[1]
		roU.SHORT[i]=roU.SHORT[i]-delta[2]		
		#print(c(iss,roL.SHORT[i],roU.SHORT[i],eq1,eq2))
		if(max(abs(delta))<eps) break
	}	
	if(curit==maxit) roL.SHORT[i]=roU.SHORT[i]=NA
}

out.typeI[4,1]=mean(roL.SHORT<ro0 & ro0<roU.SHORT,na.rm=T)	
print(Sys.time()-t0)	

out.typeI

}
