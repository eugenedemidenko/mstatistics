nrEXPcr <-
function(sigma2=0.5,b1.true=1,b2.true=.1,x=c(-4,-1,2,3,4,5),N=200,ss=24)
{
dump("nrEXPcr","c:\\Projects\\Mode\\nrEXPcr.r")
set.seed(ss)
par(mfrow=c(1,2),mar=c(4.5,4.5,1,1),cex.lab=1.5)
n=length(x)
m=2
confL=c(0.95,.8,.5)
qF=qf(confL,df1=m,df2=n-m)
b1=seq(from=-.5,to=2,length=N) #limits for beta1
b2=seq(from=-1.5,to=1.25,length=N) #limits for beta2
In=diag(rep(1,n),n,n)
for(ig in 1:2)
{
	eps=rnorm(n,mean=0,sd=sqrt(sigma2))			
	y=b1.true*exp(b2.true*x)+eps
	Q=matrix(ncol=N,nrow=N)
	for(i1 in 1:N)
	for(i2 in 1:N)
	{	
		res=y-b1[i1]*exp(b2[i2]*x)
		G=cbind(exp(b2[i2]*x),x*b1[i1]*exp(b2[i2]*x))
		iGG=solve(t(G)%*%G)
		H=G%*%iGG%*%t(G)
		Q1=t(res)%*%H%*%res
		Q2=t(res)%*%(In-H)%*%res
		Q[i1,i2]=(Q1/m)/(Q2/(n-m))
	}
	contour(b1,b2,Q,levels=qF,lty=1:3,lwd=2,drawlabels=F,xlab="beta1",ylab="beta2")
	points(b1.true,b2.true,cex=1.7)
	legend("bottomright",paste("confL =",confL),lty=1:3,cex=1.6,bg="gray95")	
}

}
