vartestOP <-
function(s20.max=2,alpha=0.05,N=100)
{ #Overall powerfull test for normal variance
dump("vartestOP","c:\\Projects\\Mode\\vartestOP.r")

#Solves the system C(q2)-C(q1)=1-alpha, (n-adj)*log(q2/q1)-(q2-q1)=0 by Newton's algorithm
var.ql=function(n,adj=1,alpha,eps=1e-05,maxit=100)
{
    n1=n-adj
    q1=qchisq(alpha/2,df=n-1)
    q2=qchisq(1-alpha/2,df=n-1)
    M=matrix(ncol=2,nrow=2)
    m=rep(NA,2)
    for(it in 1:maxit)
    {
       m[1]=n1*(log(q2)-log(q1))-(q2-q1)
       m[2]=pchisq(q1,df=n-1)-pchisq(q2,df=n-1)+(1-alpha)
       M[1,1]=n1/q1-1
       M[1,2]=-n1/q2+1
       M[2,1]=-dchisq(q1,df=n-1)
       M[2,2]=dchisq(q2,df=n-1)
       delta=solve(M)%*%m
       if(max(abs(delta))<eps)break
       q1=q1+delta[1]
       q2=q2+delta[2]
    }
	return(c(q1,q2))
}
funpow.above=function(x,s20,n,qU,qL) pchisq(qU*s20/x,df=n-1)-pchisq(qL*s20/x,df=n-1)
ffunpow.above=function(x,s20,n,qU,qL) x*dchisq(qU*s20/x,df=n-1)-x*dchisq(qL*s20/x,df=n-1)

s20=seq(from=0.1,to=s20.max,length=N)
pow.above=matrix(ncol=3,nrow=N)
par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),cex.main=1.5,cex.lab=1.5)
for(n in c(5,20))
{

	#Equal-tailed
	q1=qchisq(alpha/2,df=n-1);q2=qchisq(1-alpha/2,df=n-1)
	for(i in 1:N)
	{
		pow.above[i,1]=integrate(funpow.above,s20=s20[i],n=n,qU=q2,qL=q1,lower=0,upper=Inf)$value
		pow.above[i,1]=integrate(ffunpow.above,s20=s20[i],n=n,qU=q2,qL=q1,lower=0,upper=Inf)$value
	}

	#Unbiased test
	q=var.ql(n=n,adj=1,alpha=alpha,eps=1e-05,maxit=100)
	for(i in 1:N)
		pow.above[i,2]=integrate(funpow.above,s20=s20[i],n=n,qU=q[2],qL=q[1],lower=0,upper=Inf)$value
	
	#DL test
	q.DL=var.ql(n=n,adj=3,alpha=alpha,eps=1e-05,maxit=100)
	for(i in 1:N)
		pow.above[i,3]=integrate(funpow.above,s20=s20[i],n=n,qU=q.DL[2],qL=q.DL[1],lower=0,upper=Inf)$value
	matplot(s20,pow.above/s20,col=1,type="l",lwd=2,xlab="Null hypothesis variance",ylab="Relative area abobe power",main=paste("n =",n))

	legend("bottomright",c("Equal-tailed","Unbiased","Density level"),lty=1:3,lwd=2,cex=1.5,bg="gray97")		
}

}
