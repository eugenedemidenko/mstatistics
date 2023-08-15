linr2_5 <-
function(p=2,n=8,N=200,Nterm=500)
{
dump("linr2_5","c:\\Projects\\Mode\\linr2_5.r")
r2fder=function(x,R2P,p,n,Nterm)
{
	delta=n*R2P/(1-R2P)
	j=1:Nterm
	T1=-0.5*exp(-delta/2)*pbeta(x,p/2,(n-p-1)/2)
	T2=-0.5*exp(-delta/2)*sum(exp(-lgamma(j+1)-j*log(2)+j*log(delta))*pbeta(x,p/2+j,(n-p-1)/2))
	T3=exp(-delta/2)*sum(exp(-lgamma(j)-j*log(2)+(j-1)*log(delta))*pbeta(x,p/2+j,(n-p-1)/2))
	der=(T1+T2+T3)*n*R2P/(1-R2P)^2
	return(der)
}
par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),cex.lab=1.4,cex.main=2)
R2=seq(from=0.01,to=0.99,length=N)
R2.MCL=rep(NA,N)
n=c(10,20)
for(n in c(10,20))
{
	R2.adj=1-(1-R2)*(n-1)/(n-p-1)		
	for(i in 1:N)
		R2.MCL[i]=optimize(r2fder,x=R2[i],p=p,n=n,Nterm=Nterm,lower=0,upper=1)$minimum
	matplot(R2,cbind(R2.MCL,R2.adj),ylim=c(-.4,1),type="l",col=1,lwd=3,xlab="Observed CoD",ylab="Adjusted and MCL CoD")
	title(paste("n = ",n))
	segments(-1,0,2,0,lty=2)
	legend("topleft",c("MCL CoD","Adjusted CoD"),col=1,lty=1:2,lwd=3,cex=1.5,bg="gray96")
}
}
