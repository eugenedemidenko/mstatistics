pvalVAR <-
function(s20=0.1^2,n=10,alpha=0.05,cx=1.5,N=20,nSim=1000000)
{
#Computes 3 p-values for normal variance: equal-tail, unbiased, DL, and integral, and checks via power function 
dump("pvalVAR","c:\\projects\\Mode\\pvalVAR.r")
t0=Sys.time()
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

#computes the solutions of the equation ln(s)-s/A=B
lnsAB=function(A,B,eps=1e-10,maxit=10)
{
	s1=exp(B)
	for(it in 1:maxit)
	{
		delta=(log(s1)-s1/A-B)/(1/s1-1/A)
		if(max(abs(delta))<eps) break
		s1=s1-delta
	}
	s2=A*(log(A)-B)
	for(it in 1:maxit)
	{
		delta=(log(s2)-s2/A-B)/(1/s2-1/A)
		if(max(abs(delta))<eps) break
		s2=s2-delta
	}
	return(cbind(s1,s2))
}
par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
s2s=seq(from=s20/3,to=s20*2,length=N)
qL.ET=qchisq(alpha/2,df=n-1)
qU.ET=qchisq(1-alpha/2,df=n-1)
pow.ET=1+pchisq(s20*qL.ET/s2s,df=n-1)-pchisq(s20*qU.ET/s2s,df=n-1)
qLU.DL=var.ql(n=n,adj=3,alpha=alpha)
pow.DL=1+pchisq(s20*qLU.DL[1]/s2s,df=n-1)-pchisq(s20*qLU.DL[2]/s2s,df=n-1)
qLU.UNB=var.ql(n=n,adj=1,alpha=alpha)
pow.UNB=1+pchisq(s20*qLU.UNB[1]/s2s,df=n-1)-pchisq(s20*qLU.UNB[2]/s2s,df=n-1)
#plotting 3 theoretical power functions
matplot(s2s,cbind(pow.ET,pow.DL,pow.UNB),col=2:4,lwd=2,lty=1:3,type="l",xlab="Alternative variance",ylab="Power, probability")
segments(s20,-1,s20,.2,lty=2)
segments(-1,alpha,100,alpha,lty=2)
legend("topright",c("ET","DL","UNB"),lty=1:3,lwd=2,pch=1,col=2:4,cex=1.5,bg="grey95")
text(s20,0.4,paste("n =",n,"\nnull variance =",s20),font=3,cex=1.5)
S0=rchisq(nSim,df=n-1)
A=(n-3)*s20
pUNB=rep(NA,nSim)
#checking theoretical powers via p-values as proportion p <= alpha
for(i in 1:N)
{
	S=S0*s2s[i]
	#Equal-tail test
	pET=rep(NA,nSim)
	p1ET=2*pchisq(S/s20,df=n-1);p2ET=2*pchisq(S/s20,df=n-1,lower.tail=F)
	pET[p1ET<=p2ET]=p1ET[p1ET<=p2ET]
	pET[p1ET>p2ET]=p2ET[p1ET>p2ET]
	points(s2s[i],mean(pET<=alpha),col=2,cex=cx)
	#Maximum concentration statistics/density level test
	qLU=var.ql(n=n,adj=3,alpha=alpha)
	div=pchisq(qLU[1],df=n-1)/alpha
	pS=pchisq(S/s20,df=n-1)
	L1=1/div*pS
	L2=alpha/(alpha-pchisq(qLU[1],df=n-1))*(1-pS)
	pDL=L1
	pDL[pS>div]=L2[pS>div]
	points(s2s[i],mean(pDL<=alpha),col=3,cex=cx)			
	#alternative via integral
	B=log(S)-S/A
	s12=lnsAB(A=A,B=B)
	powDL_INT=1-pchisq(s12[,2]/s20,df=n-1)+pchisq(s12[,1]/s20,df=n-1)
	points(s2s[i],mean(powDL_INT<=alpha),col=3,pch=2,cex=cx)			
	#Unbiased test
	qLU=var.ql(n=n,adj=1,alpha=alpha)
	div=pchisq(qLU[1],df=n-1)/alpha
	L1=1/div*pS
	L2=alpha/(alpha-pchisq(qLU[1],df=n-1))*(1-pS)
	pUNB=L1
	pUNB[pS>div]=L2[pS>div]
	points(s2s[i],mean(pUNB<=alpha),col=4,cex=cx)					
	print(c(mean(pDL<=alpha),mean(powDL_INT<=alpha)))
}
Sys.time()-t0

}
