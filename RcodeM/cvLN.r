cvLN <-
function(n=10,alpha=0.05,ksi.true=.3,nSim=100000,N=100)
{
dump("cvLN","c:\\Projects\\Mode\\cvLN.r")
qLU=var.ql(n=n,adj=1,alpha=alpha)
ksi.alt=seq(from=.1,to=.75,length=N)
ln=log(1+ksi.true^2)/log(1+ksi.alt^2)
Pksi=1+pchisq(qLU[1]*ln,df=n-1)-pchisq(qLU[2]*ln,df=n-1)
par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
plot(ksi.alt,Pksi,lwd=2,type="l",ylim=c(0,1),xlab="Alternative coefficient of variation",ylab="Power, probability")
segments(-1,alpha,1000,alpha)
segments(ksi.true,-1,ksi.true,1)
text(.4,.8,paste("n =",n),cex=1.75,font=3)
Pksi0=1+pchisq(qLU[1],df=n-1)-pchisq(qLU[2],df=n-1)
points(ksi.true,Pksi0,cex=1.5)
iN=seq(from=1,to=N,length=N/2)
Z1a=qnorm(1-alpha/2)
as.pow=rep(NA,N/2)
for(i in 1:(N/2))
{
	sigma2=log(1+ksi.alt[iN[i]]^2)
	S=sigma2*rchisq(nSim,df=n-1)
	s2L=S/qLU[2];s2U=S/qLU[1]
	ksiL=sqrt(exp(s2L)-1);ksiU=sqrt(exp(s2U)-1)
	pow=mean(ksiL>ksi.true | ksiU<ksi.true)
	points(ksi.alt[iN[i]],pow,cex=1.25)
	
	s2.est=S/(n-1)
	es2=exp(s2.est)
	ksi.est=sqrt(es2-1)
	var.ksi=s2.est^2*es2^2/2/(n-1)/(es2-1)
	as.L=ksi.est-Z1a*sqrt(var.ksi)
	as.U=ksi.est+Z1a*sqrt(var.ksi)
	as.pow[i]=mean(as.L>ksi.true | as.U<ksi.true)
	points(ksi.alt[iN[i]],as.pow[i],pch=2)
}
legend("topright",c("Exact unbiased","Asymptotic"),lty=c(1,0),lwd=c(2,NULL),pch=c(1,2),cex=1.6,bg="gray97")

}
