sdM <-
function(n=10,sigma.true=1,s.CI=.85,s.HT=1.2,s.PV=1.1,alpha=0.05,N=200,nSim=1000000,dr="c")
{
dump("sdM",paste(dr,":\\Projects\\Mode\\sdM.r",sep=""))


SR=rchisq(nSim,df=n-1)
SR.CI=s.CI^2*SR #S under the alternative s.CI
SR.HT=s.HT^2*SR #S under the alternative s.HT
SR.PV=s.PV^2*SR

par(mfrow=c(1,1),mar=c(4,4,1,1))
sds=seq(from=.7,to=1.3,length=N)
symb=c(1,2,5,7)

#Equal-tailed inference
q1=qchisq(alpha/2,df=n-1);q2=qchisq(1-alpha/2,df=n-1)
pow.eqtail=1+pchisq(q1*sigma.true^2/sds^2,df=n-1)-pchisq(q2*sigma.true^2/sds^2,df=n-1)
plot(sds,pow.eqtail,type="l",lwd=3,ylab="",xlab="",ylim=c(0.03,.3))
segments(s.CI,-1,s.CI,.14,col=grey(.8));text(s.CI,.141,"Confidence interval",srt=90,adj=0,cex=1.25,font=3)
segments(s.HT,-1,s.HT,.23,col=grey(.8));text(s.HT,.232,"Rejection rule",srt=90,adj=0,cex=1.25,font=3)
segments(s.PV,-1,s.PV,.14,col=grey(.8));text(s.PV,.141,"P-value rejection",srt=90,adj=0,cex=1.25,font=3)
mtext(side=1,"Alternative standard deviation",cex=1.5,line=2.75)
mtext(side=2,"Power, probability",cex=1.5,line=2.5)
segments(-1,.05,3,.05);segments(sigma.true,-1,sigma.true,.08)
pow.eqtail.sim.HT=1-mean(sigma.true<sqrt(SR.HT/q1) & sqrt(SR.HT/q2)<sigma.true)
points(s.HT,pow.eqtail.sim.HT,cex=2,pch=symb[1])
pow.eqtail.sim.CI=1-mean(sqrt(SR.CI/q1)>sigma.true & sqrt(SR.CI/q2)<sigma.true)
points(s.CI,pow.eqtail.sim.CI,cex=2,pch=symb[1])
pv=2*pmin(pchisq(SR.PV/sigma.true^2,df=n-1),pchisq(SR.PV/sigma.true^2,df=n-1,lower.tail=F))
points(s.PV,mean(pv<=alpha),cex=2,pch=symb[1])

#MO statistics: unbiased inference
q=var.ql(n=n,adj=1,alpha=alpha,eps=1e-05,maxit=100)
pow.dl=1+pchisq(q[1]*sigma.true^2/sds^2,df=n-1)-pchisq(q[2]*sigma.true^2/sds^2,df=n-1)
lines(sds,pow.dl,lwd=3,lty=2)
pow.ql.sim.HT=1-mean(sigma.true<sqrt(SR.HT/q[1]) & sqrt(SR.HT/q[2])<sigma.true)
points(s.HT,pow.ql.sim.HT,cex=2,pch=symb[2])
pow.ql.sim.CI=1-mean(sigma.true<sqrt(SR.CI/q[1]) & sqrt(SR.CI/q[2])<sigma.true)
points(s.CI,pow.ql.sim.CI,cex=2,pch=symb[2])

pS=pchisq(SR.PV/sigma.true^2,df=n-1)
div=pchisq(q[1],df=n-1)/alpha
L1=1/div*pS
L2=alpha/(alpha-pchisq(q[1],df=n-1))*(1-pS)
pUNB=L1
pUNB[pS>div]=L2[pS>div]
points(s.PV,mean(pUNB<=alpha),cex=2,pch=symb[2])


#MC statistics: short inference
q=var.ql(n=n,adj=3,alpha=alpha,eps=1e-05,maxit=100)
pow.dl=1+pchisq(q[1]*sigma.true^2/sds^2,df=n-1)-pchisq(q[2]*sigma.true^2/sds^2,df=n-1)
lines(sds,pow.dl,lwd=3,lty=3)
pow.ql.sim.HT=1-mean(sigma.true<sqrt(SR.HT/q[1]) & sqrt(SR.HT/q[2])<sigma.true)
points(s.HT,pow.ql.sim.HT,cex=2,pch=symb[3])
pow.ql.sim.CI=1-mean(sigma.true<sqrt(SR.CI/q[1]) & sqrt(SR.CI/q[2])<sigma.true)
points(s.CI,pow.ql.sim.CI,cex=2,pch=symb[3])

div=pchisq(q[1],df=n-1)/alpha
L1=1/div*pS
L2=alpha/(alpha-pchisq(q[1],df=n-1))*(1-pS)
pDL=L1
pDL[pS>div]=L2[pS>div]
points(s.PV,mean(pDL<=alpha),cex=2,pch=symb[3])

#Test dual to short CI
q=var.ql(n=n,adj=0,alpha=alpha,eps=1e-05,maxit=100)
pow.dl=1+pchisq(q[1]*sigma.true^2/sds^2,df=n-1)-pchisq(q[2]*sigma.true^2/sds^2,df=n-1)
lines(sds,pow.dl,lwd=3,lty=3)
pow.ql.sim.HT=1-mean(sigma.true<sqrt(SR.HT/q[1]) & sqrt(SR.HT/q[2])<sigma.true)
points(s.HT,pow.ql.sim.HT,cex=2,pch=symb[4])
pow.ql.sim.CI=1-mean(sigma.true<sqrt(SR.CI/q[1]) & sqrt(SR.CI/q[2])<sigma.true)
points(s.CI,pow.ql.sim.CI,cex=2,pch=symb[4])

div=pchisq(q[1],df=n-1)/alpha
L1=1/div*pS
L2=alpha/(alpha-pchisq(q[1],df=n-1))*(1-pS)
pDL=L1
pDL[pS>div]=L2[pS>div]
points(s.PV,mean(pDL<=alpha),cex=2,pch=symb[4])


legend("topleft",c("Equal-tail (traditional)","Unbiased (MO statistic)","Density level (MC statistic)","Test dual to short CI"),lty=1:4,pch=symb,lwd=3,cex=1.8,bg="gray97")
text(sigma.true,.13,paste("n =",n),cex=2,font=3)
}
