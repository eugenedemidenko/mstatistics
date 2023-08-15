pvalVAR_comp <-
function(job=1,n=10,alpha=0.05,N=2000)
{
#Computes 3 p-values for normal variance: equal-tail, unbiased, DL, and integral, and checks via power function 
dump("pvalVAR_comp","c:\\projects\\Mode\\pvalVAR_comp.r")

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
if(job==1)
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	S=seq(from=0.001,to=2.5*n,length=N)
	alphas=c(0.01,0.05,0.2,.5);nal=length(alphas)
	p_alph=matrix(ncol=nal,nrow=N)
	for(i in 1:nal)
	{
		qL=var.ql(n=n,adj=3,alpha=alphas[i])[1]
		div=pchisq(qL,df=n-1)/alphas[i]
		pS=pchisq(S,df=n-1)
		L1=pS/div
		L2=alphas[i]/(alphas[i]-pchisq(qL,df=n-1))*(1-pS)
		pDL=L1
		pDL[pS>div]=L2[pS>div]
		p_alph[,i]=pDL			
	}
	A=n-3
	B=log(S)-S/A
	s12=lnsAB(A=A,B=B)
	pDL_alt=1-abs(pchisq(s12[,2],df=n-1)-pchisq(s12[,1],df=n-1))
	matplot(pDL_alt,p_alph,type="l",ylim=c(0,1),col=1,lwd=2,xlab="Alternative p-value",ylab="Alpha-dependent p-value")
	text(0.05,.9,paste("n =",n),cex=1.5,adj=0,font=3)
	for(i in 1:nal) 
	text(pDL_alt[which(p_alph[,i]==max(p_alph[,i]))],rep(.96,nal)+.06,as.character(i),font=2)
	legend("bottomright",paste(1:nal,": alpha=",alphas,sep=""),lty=1:nal,lwd=2,cex=1.5,bg="grey96")
}
if(job==2)
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	S=pET=seq(from=0.001,to=2*n,length=N)
#ET p-value
	p1ET=2*pchisq(S,df=n-1);p2ET=2*pchisq(S,df=n-1,lower.tail=F)
	pET[p1ET<=p2ET]=p1ET[p1ET<=p2ET]
	pET[p1ET>p2ET]=p2ET[p1ET>p2ET]
	etaET=qchisq(0.5,df=n-1)/(n-1)
	
#alpha-dependent DL p-value  
	qLU=var.ql(n=n,adj=3,alpha=alpha) #DL p-value
	div=pchisq(qLU[1],df=n-1)/alpha
	pS=pchisq(S,df=n-1)
	L1=1/div*pS
	L2=alpha/(alpha-pchisq(qLU[1],df=n-1))*(1-pS)
	pDL=L1
	pDL[pS>div]=L2[pS>div]
	etaDL=1/(n-1)*qchisq(1/alpha*pchisq(qLU[1],df=n-1),df=n-1)
#alternative DL p-value
	A=(n-3)
	B=log(S)-S/A
	s12=lnsAB(A=A,B=B)
	pDL_alt=1-abs(pchisq(s12[,2],df=n-1)-pchisq(s12[,1],df=n-1))
	etaADL=(n-3)/(n-1)
#Unbiased p-value
	qLU=var.ql(n=n,adj=1,alpha=alpha)
	div=pchisq(qLU[1],df=n-1)/alpha
	L1=1/div*pS
	L2=alpha/(alpha-pchisq(qLU[1],df=n-1))*(1-pS)
	pUNB=L1
	pUNB[pS>div]=L2[pS>div]
	
	etaUNB=1/(n-1)*qchisq(1/alpha*pchisq(qLU[1],df=n-1),df=n-1)
	
	allp=cbind(pET,pDL_alt,pDL,pUNB)
	matplot(S/(n-1),allp,type="l",col=1,lwd=2,ylab="P-value",xlab="Normalized variance")
	segments(-1,alpha,max(S),alpha,lty=2)
	text(1.35,alpha+.03,paste("a =",alpha),font=5,cex=1.5)
	text(0.01,0.8,paste("n =",n),cex=1.5,font=3,adj=0)
	segments(etaET,-1,etaET,1,lty=2)
	segments(etaDL,-1,etaDL,1,lty=2)
	segments(etaADL,-1,etaADL,1,lty=2)
	segments(etaUNB,-1,etaUNB,1,lty=2)
	
	for(i in 1:4) 
	text((S/(n-1))[which(allp[,i]==max(allp[,i]))],rep(.96,4)+.06,as.character(i),font=2)
	legend("topright",c("1: ET","2: ADL","3: DL","4: UNB"),lty=1:4,lwd=2,cex=1.25,bg="grey95")
}
}
