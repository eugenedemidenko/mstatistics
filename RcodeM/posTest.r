posTest <-
function(job=1,lambda0=1.5,alpha=0.05,N=200,nSim=1000000)
{# hypothesis testing for Poisson rate parameter lambda
dump("posTest","c:\\Projects\\Mode\\posTest.r")
#install.packages("nleqslv")
library(nleqslv)
sys2=function(qLU,n,lambda,alpha)
{
	qL=qLU[1];qU=qLU[2]
	eq1=pchisq(2*lambda*n,df=2*(qL+1))-pchisq(2*lambda*n,df=2*(qU+1))-(1-alpha)
	eq2=dchisq(2*lambda*n,df=2*(qL+1))-dchisq(2*lambda*n,df=2*(qU+1))
	return(c(eq1,eq2))
}
if(job==1) #Equal-tail and unbiased tests
{
	par(mfrow=c(1,2),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	for(n in c(2,5))
	{
		lambdas=seq(from=0,to=5*lambda0,length=N)
		ql_ET=qpois(alpha/2,lambda=n*lambda0)	
		qu_ET=qpois(1-alpha/2,lambda=n*lambda0)	
		q_UNB=nleqslv(c(ql_ET,qu_ET),sys2,n=n,lambda=lambda0,alpha=alpha)$x
		ql_UNB=max(0,q_UNB[1]);qu_UNB=q_UNB[2]

		powET=1+ppois(ql_ET,lambda=n*lambdas)-ppois(qu_ET,lambda=n*lambdas)
		powUNB=1+ppois(ql_UNB,lambda=n*lambdas)-ppois(qu_UNB,lambda=n*lambdas)
				
		matplot(lambdas,cbind(powET,powUNB),type="l",col=1,lwd=2,ylim=c(0,1),xlab="Alternative lambda",ylab="Power function, probability")
		segments(-1,alpha,100,alpha,lty=2)
		segments(lambda0,-1,lambda0,1,lty=2)
		text(6,.5,paste("l =",lambda0),font=5,cex=2)
		text(6,.4,paste("n =",n),font=3,cex=2)
		for(i in 1:(N/10))
		{
			lambda=lambdas[i*10]
			S=rpois(nSim,lambda=lambda*n)
			powSIM_ET=1-mean(S<=qu_ET & S>ql_ET)
			points(lambda,powSIM_ET,cex=1.5)
			powSIM_UNB=1-mean(S<=qu_UNB & S>ql_UNB)
			points(lambda,powSIM_UNB,pch=2,cex=1.5)
		}
		legend("bottomright",c("Equal-tail","Unbiased"),lty=1:2,lwd=2,pch=1:2,cex=1.75,bg="gray96")
	}
}
if(job==2)# computation of the p-value for the example
{
	lambda0=1.2
	S=11;n=5
	F_left=ppois(S,lambda=n*lambda0)	
	F_right=ppois(S,lambda=n*lambda0,lower.tail=F)	
	p_ET=2*min(F_left,F_right)
	cat("Equal-tail test p-value:",round(p_ET,5))
	
	ql_ET=qpois(alpha/2,lambda=n*lambda0)	
	qu_ET=qpois(1-alpha/2,lambda=n*lambda0)	
	cat("\nQuantiles of the equal-tail test:",ql_ET,qu_ET)
	q_UNB=nleqslv(c(ql_ET,qu_ET),sys2,n=n,lambda=lambda0,alpha=alpha)$x
	cat("\nQuantiles of the unbiased test:",round(q_UNB[1],5),round(q_UNB[2],5))
	FS0=ppois(S,lambda=n*lambda0)	
	FS0_r=ppois(S,lambda=n*lambda0,lower.tail=F)	
	FL0=ppois(q_UNB[1],lambda=n*lambda0)	
	if(FS0<=FL0/alpha) pv_UNB=alpha/FL0*FS0 else pv_UNB=alpha/(alpha-FL0)*FS0_r
	cat("\nUnbiased test p-value:",round(pv_UNB,5),"\n")
	
}

}
