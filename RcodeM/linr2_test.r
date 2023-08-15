linr2_test <-
function(job=1,s2=.1,m=2,n=8,r2.est=.4,N=200,Nterm=50,lambda=0.95)
{
dump("linr2_test","c:\\Projects\\Mode\\linr2_test.r")
#install.packages("hypergeo")
print(date())
library(hypergeo)
#install.packages("nleqslv")
library(nleqslv)
fder=function(x,R2P,m,n,Nterm)
{
	delta=n*R2P/(1-R2P)
	j=1:Nterm
	T1=-0.5*exp(-delta/2)*pbeta(x,m/2,(n-m-1)/2)
	T2=-0.5*exp(-delta/2)*sum(exp(-lgamma(j+1)-j*log(2)+j*log(delta))*pbeta(x,m/2+j,(n-m-1)/2))
	T3=exp(-delta/2)*sum(exp(-lgamma(j)-j*log(2)+(j-1)*log(delta))*pbeta(x,m/2+j,(n-m-1)/2))
	der=(T1+T2+T3)*n/(1-R2P)^2
	return(der)
}
eqtest=function(RLU,m,n,R2P0,Nterm,lambda=0.95)
{
	delta=n*R2P0/(1-R2P0)
	RL=RLU[1];RU=RLU[2]
	FL=pf(RL/(1-RL)*(n-m-1)/m,df1=m,df2=n-m-1,ncp=delta)
	FU=pf(RU/(1-RU)*(n-m-1)/m,df1=m,df2=n-m-1,ncp=delta)
	F1L=fder(x=RL,R2P=R2P0,m=m,n=n,Nterm=Nterm)
	F1U=fder(x=RU,R2P=R2P0,m=m,n=n,Nterm=Nterm)
	return(c(FU-FL-lambda,F1U-F1L))		
}


R2P0=r2.est
delta=n*R2P0/(1-R2P0)
q=qf((1-lambda)/2,df1=m,df2=n-m-1,ncp=delta)
RL0=q*m/(n-m-1)/(1+q*m/(n-m-1))
q=qf((1+lambda)/2,df1=m,df2=n-m-1,ncp=delta)
RU0=q*m/(n-m-1)/(1+q*m/(n-m-1))
RLU=c(RL0,RU0)
#print(RLU)
out=nleqslv(RLU,eqtest,m=m,n=n,R2P0=R2P0,Nterm=Nterm,lambda=lambda)
#print(out)
RL=out$x[1];RU=out$x[2]
par(mfrow=c(1,1))
R2=seq(from=0,to=1,length=N)
delta=n*R2/(1-R2)
Ppr=1+pf(RL/(1-RL)*(n-m-1)/m,df1=m,df2=n-m-1,ncp=delta)-pf(RU/(1-RU)*(n-m-1)/m,df1=m,df2=n-m-1,ncp=delta)
par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
plot(R2,Ppr,type="l",lwd=3,xlab="Alternative parent CoD",ylab="Power, probability")
segments(-1,1-lambda,1,1-lambda,col=2)
segments(R2P0,-1,R2P0,.4,col=3)

Ppr.eq=1+pf(RL0/(1-RL0)*(n-m-1)/m,df1=m,df2=n-m-1,ncp=delta)-pf(RU0/(1-RU0)*(n-m-1)/m,df1=m,df2=n-m-1,ncp=delta)
lines(R2,Ppr.eq)	
legend("topleft",c("Unbiased test","Equal-tail biased test"),lwd=c(3,1),lty=1,cex=1.5,bg="gray95")
text(0.9,(1-lambda)+.04,"a = 0.05",font=5,cex=1.5)	

}
