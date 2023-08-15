qual_cont <-
function(s20=0.1,alpha=0.05)
{
dump("qual_cont","c:\\projects\\Mode\\qual_cont.r")
library(EnvStats)

#computes qL for DL and UNB tests
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
lnsAB=function(A,B,eps=1e-05,maxit=100)
{
	s1=exp(B)
	for(it in 1:maxit)
	{
		delta=(log(s1)-s1/A-B)/(1/s1-1/A)
		if(abs(delta)<eps) break
		s1=s1-delta
	}
	
	s2=A*(log(A)-B)
	for(it in 1:maxit)
	{
		delta=(log(s2)-s2/A-B)/(1/s2-1/A)
		if(abs(delta)<eps) break
		s2=s2-delta
	}
	return(c(s1,s2))
}

Y=c(7.34, 6.95, 7.73, 7.23, 6.71, 7.85)
n=length(Y)
pvals=as.data.frame(rep(NA,5),row.names=c("EnvStats","Equal-tail","DL-alpha","DL-alternative","Unbiased"))
names(pvals)="p-value"
#varTest function of of EnvStats
p.varTest=varTest(x=Y,sigma.squared=s20)$p.value
pvals[1,1]=p.varTest

S=var(Y)*(n-1)
#Equal-tail
cdfL=pchisq(S/s20,df=n-1)
cdfU=pchisq(S/s20,df=n-1,lower.tail=F)
p.ET=2*min(cdfL,cdfU)
pvals[2,1]=p.ET

#DL-alpha
qL=var.ql(n=n,adj=3,alpha=alpha)[1]
div=pchisq(qL,df=n-1)/alpha
pS=pchisq(S/s20,df=n-1)
L1=pS/div
L2=alpha/(alpha-pchisq(qL,df=n-1))*(1-pS)
pDL=L1
pDL[pS>div]=L2[pS>div]
pvals[3,1]=pDL

#DL-alternative
A=(n-3)*s20
B=log(S)-S/A
s12=lnsAB(A=A,B=B)
pDL_alt=1-abs(pchisq(s12[1]/s20,df=n-1)-pchisq(s12[2]/s20,df=n-1))
pvals[4,1]=pDL_alt

#Unbiased
qL=var.ql(n=n,adj=1,alpha=alpha)[1]
div=pchisq(qL,df=n-1)/alpha
pS=pchisq(S/s20,df=n-1)
L1=1/div*pS
L2=alpha/(alpha-pchisq(qL,df=n-1))*(1-pS)
pUNB=L1
pUNB[pS>div]=L2[pS>div]
pvals[5,1]=pUNB

pvals
}
