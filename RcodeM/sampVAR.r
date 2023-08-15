sampVAR <-
function(var.0=2,var.alt=3,p=0.8,alpha=0.05)
{
	dump("sampVAR","c:\\Projects\\Mode\\sampVAR.r")
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

	n.apr=(qnorm(1-alpha/2)+qnorm(p))^2*2*var.alt^2/(var.alt-var.0)^2
	n=(round(n.apr/3)):round((n.apr*3))
	L=length(n)
	pow=rep(NA,L)
	for(i in 1:L)
	{
		qLU=var.ql(n=n[i],adj=1,alpha=alpha)
		pow[i]=1+pchisq(qLU[1]*var.0/var.alt,df=n[i]-1)-pchisq(qLU[2]*var.0/var.alt,df=n[i]-1)	
	}
	par(mfrow=c(1,1),mar=c(4.5,4.5,3,1),cex.lab=1.5,cex.main=1.5)
	plot(n,pow,type="l",lwd=3,xlab="Sample size, n",ylab="Power",main=paste("Required unbiased power, p =",p,", var.0 =",var.0,", var.alt =",var.alt))
	segments(0,p,max(n),p,lty=2)
	segments(n.apr,-1,n.apr,2,lty=2)
	a=abs(pow-p)
	n.req=n[a==min(a)]
	segments(n.req,-1,n.req,2)
	legend("bottomright",c(paste("n appr =",round(n.apr)),paste("n exact =",round(n.req))),col=1,lty=c(2,1),pch=c(1,16),lwd=2,cex=1.75,bg="gray95")
	
	Pn=pnorm(sqrt(n)*(var.alt-var.0)/sqrt(2)/var.alt-qnorm(1-alpha/2))
	lines(n,Pn,lty=2,lwd=2)
	Pnit=pnorm(sqrt(n.apr)*(var.alt-var.0)/sqrt(2)/var.alt-qnorm(1-alpha/2))
	points(n.apr,Pnit,pch=1,cex=1.5)
#Modified Newton's algotithm	
	nit=n.apr
	rel=(var.alt-var.0)/sqrt(2)/var.alt
	derp.apr=dnorm(sqrt(nit)*rel-sign(var.alt-var.0)*qnorm(1-alpha/2))*rel/2/sqrt(nit)
	for(it in 1:10)
	{
		qLU=var.ql(n=nit,adj=1,alpha=alpha)
		pnit=1+pchisq(qLU[1]*var.0/var.alt,df=nit-1)-pchisq(qLU[2]*var.0/var.alt,df=nit-1)
		delta=(pnit-p)/derp.apr*sign(var.alt-var.0)	
		nit=nit-delta
		if(abs(delta)<0.01) break
		#print(c(it,nit,pnit))
	}
	points(nit,pnit,pch=16,cex=1.5)
	cat("n exact =",round(nit),"\n")

}
