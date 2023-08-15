vartest <-
function(n=20,s20=1,s2.CI=.85,s2.HT=1.2,alpha=0.05,N=200,nSim=10000000,dr="c")
{ #power function and CI for normal variance
dump("vartest",paste(dr,":\\Projects\\Mode\\vartest.r",sep=""))
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

SR=rchisq(nSim,df=n-1)
SR.CI=s2.CI*SR #S under the alternative s2.CI
SR.HT=s2.HT*SR #S under the alternative s2.HT

symb=c(1,2,5,7)
par(mfrow=c(1,1),mar=c(4,4,1,1))
s2=seq(from=.7,to=1.5,length=N)

#Equal-tailed inference
q1=qchisq(alpha/2,df=n-1);q2=qchisq(1-alpha/2,df=n-1)
pow.eqtail=1+pchisq(q1*s20/s2,df=n-1)-pchisq(q2*s20/s2,df=n-1) #theoretical power
plot(s2,pow.eqtail,type="l",lwd=3,ylab="",xlab="",ylim=c(0,.3))
segments(s2.CI,-1,s2.CI,.12,col=grey(.8));text(s2.CI,.121,"Confidence interval",srt=90,adj=0,cex=1.25,font=3)
segments(s2.HT,-1,s2.HT,.15,col=grey(.8));text(s2.HT,.151,"Rejection rule",srt=90,adj=0,cex=1.25,font=3)
mtext(side=1,"Alternative variance",cex=1.5,line=2.5)
mtext(side=2,"Power, probability",cex=1.5,line=2.75)
segments(-1,.05,3,.05);segments(s20,-1,s20,.08);text(s20,.081,"Null variance",cex=1.25,srt=90,font=3,adj=0)
pow.eqtail.sim.HT=1-mean(q1<SR.HT/s20 & SR.HT/s20<q2) #simulation-derived power from rejection rule
points(s2.HT,pow.eqtail.sim.HT,cex=2,pch=symb[1])
pow.eqtail.sim.CI=1-mean(SR.CI/q1>s20 & SR.CI/q2<s20) #simulation-derived power from rejection by CI
points(s2.CI,pow.eqtail.sim.CI,cex=2,pch=symb[1])

#MO statistics: unbiased inference
q=var.ql(n=n,adj=1,alpha=alpha,eps=1e-05,maxit=100)
pow.unb=1+pchisq(q[1]*s20/s2,df=n-1)-pchisq(q[2]*s20/s2,df=n-1) #theoretical power
lines(s2,pow.unb,lwd=3,lty=2)
pow.unb.sim.HT=1-mean(q[1]<SR.HT/s20 & SR.HT/s20<q[2]) #simulation-derived power from rejection rule
points(s2.HT,pow.unb.sim.HT,cex=2,pch=symb[2])
pow.unb.sim.CI=1-mean(SR.CI/q[1]>s20 & SR.CI/q[2]<s20) #simulation-derived power from rejection by CI
points(s2.CI,pow.unb.sim.CI,cex=2,pch=symb[2])

#MC statistics: DL test
q.DL=var.ql(n=n,adj=3,alpha=alpha,eps=1e-05,maxit=100)
pow.dl=1+pchisq(q.DL[1]*s20/s2,df=n-1)-pchisq(q.DL[2]*s20/s2,df=n-1) #theoretical power
lines(s2,pow.dl,lwd=3,lty=3)
pow.ql.sim.HT=1-mean(q.DL[1]<SR.HT/s20 & SR.HT/s20<q.DL[2]) #simulation-derived power from rejection rule
points(s2.HT,pow.ql.sim.HT,cex=2,pch=symb[3])
q.CI=var.ql(n=n,adj=3,alpha=alpha,eps=1e-05,maxit=100)
pow.ql.sim.CI=1-mean(SR.CI/q.CI[1]>s20 & SR.CI/q.CI[2]<s20) #simulation-derived power from rejection by CI
points(s2.CI,pow.ql.sim.CI,cex=2,pch=symb[3])
s2.hat=SR.HT/(n-1) 


#MC statistics: Test dual to short CI
q.DL=var.ql(n=n,adj=-1,alpha=alpha,eps=1e-05,maxit=100)
pow.dl=1+pchisq(q.DL[1]*s20/s2,df=n-1)-pchisq(q.DL[2]*s20/s2,df=n-1) #theoretical power
lines(s2,pow.dl,lwd=3,lty=3)
pow.ql.sim.HT=1-mean(q.DL[1]<SR.HT/s20 & SR.HT/s20<q.DL[2]) #simulation-derived power from rejection rule
points(s2.HT,pow.ql.sim.HT,cex=2,pch=symb[4])
q.CI=var.ql(n=n,adj=-1,alpha=alpha,eps=1e-05,maxit=100)
pow.ql.sim.CI=1-mean(SR.CI/q.CI[1]>s20 & SR.CI/q.CI[2]<s20) #simulation-derived power from rejection by CI
points(s2.CI,pow.ql.sim.CI,cex=2,pch=symb[4])
s2.hat=SR.HT/(n-1) 
legend("topleft",c("Equal-tailed (traditional)","Unbiased (MO statistic)","Density level (MC statistic)","Test dual to short CI"),lty=1:4,pch=symb,lwd=3,cex=1.8,bg="gray97")


}
