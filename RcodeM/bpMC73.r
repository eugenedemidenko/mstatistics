bpMC73 <-
function(job=1,lambda=0.95,N=100)
{
dump("bpMC73","c:\\Projects\\Mode\\bpMC73.r")
alpha=1-lambda
if(job==1)
{
	n1=10;p1=0.2;n2=20;p2=0.6
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	m1=0:n1;m2=0:n2
	P=matrix(nrow=n1+1,ncol=n2+1)
	for(i in 0:n1) 
		P[i+1,]=dbinom(i,size=n1,prob=p1)*dbinom(m2,size=n2,prob=p2)
	Pmax=max(P)
	kap0=Pmax*alpha
	print(paste("kappa apprx =",round(kap0,4)))
	kaps=unique(as.vector(P))
	kaps=sort(kaps);nk=length(kaps)
	sumk=rep(NA,nk)
	for(i in 1:nk)
		sumk[i]=sum(P[P>=kaps[i]])
	plot(kaps,sumk,lwd=2,type="s",xlab="",ylab="Probability",main="")
	mtext(side=1,"k",font=5,line=2.5,cex=2)
	segments(-1,lambda,100,lambda)
	kap=max(kaps[sumk >= lambda])
	segments(kap,sumk[kaps==kap],kap,-1)	
	segments(kap0,max(sumk[kaps >= kap0]),kap0,-1,lty=2)
	text(.01,.2,paste("k =",round(kap,4)),cex=1.75,font=5,adj=0)
	text(.02,.98,paste("l =",lambda),font=5,cex=1.5,adj=0)
	legend("topright",c("exact","approximate"),lty=1:2,cex=1.5,bg="gray96")
}
if(job==2)
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	n1=10;p1=0.2;n2=20;p2=0.6
	m1=0:n1;m2=0:n2
	P=matrix(nrow=n1+1,ncol=n2+1)
	for(i in 0:n1) 
		P[i+1,]=dbinom(i,size=n1,prob=p1)*dbinom(m2,size=n2,prob=p2)
	kap=0.0037
	plot(1,1,xlim=c(0,n1),ylim=c(0,n2),xlab="i",ylab="j",type="n")
	for(i in 0:n1)
	for(j in 0:n2)
		if(P[i+1,j+1]>=kap) points(i,j,pch=16)	
	m1=5;m2=9
	points(m1,m2,cex=1.5,pch=2)
	p12=dbinom(m1,size=n1,prob=p1)*dbinom(m2,size=n2,prob=p2)
	pval=sum(P[P<p12])
	text(5.7,9,paste("p-value =",round(pval,3)),cex=1.25,font=3,adj=0)	
	qchi=qchisq(lambda,df=2)
	for(i in 0:n1)
	for(j in 0:n2)
	{
		lhs=n1*(i/n1-p1)^2/p1/(1-p1)+n2*(j/n2-p2)^2/p2/(1-p2)
		if(lhs<qchi) points(i+.1,j+.1,col=2)	
	}
	legend("bottomright",c("Exact DL test","Wald test"),pch=c(16,1),col=c(1,2),cex=1.5,bg="gray95")
}
if(job==3)
{
	t0=Sys.time()
	n1=10;n2=20;m1.obs=2;m2.obs=4
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)	
	p1=p2=seq(from=0,to=1,length=N)
	kap=Pm=matrix(ncol=N,nrow=N)
	P=matrix(nrow=n1+1,ncol=n2+1)
	m2=0:n2	
	for(ii in 1:N)
	for(jj in 1:N)
	{
		for(i in 0:n1) 
			P[i+1,]=dbinom(i,size=n1,prob=p1[ii])*dbinom(m2,size=n2,prob=p2[jj])
		Pmax=max(P)
		kaps=unique(as.vector(P))
		kaps=sort(kaps);nk=length(kaps)
		sumk=rep(NA,nk)
		for(i in 1:nk)
			sumk[i]=sum(P[P>=kaps[i]])
		kap[ii,jj]=max(kaps[sumk >= lambda])				
		Pm[ii,jj]=dbinom(m1.obs,size=n1,prob=p1[ii])*dbinom(m2.obs,size=n2,prob=p2[jj])
	}	
	contour(p1,p2,Pm-kap,lwd=2,levels=0,xlim=c(0,.7),ylim=c(0,.7),xlab="p1",ylab="p2",drawlabels=F)
	p12=p1%*%t(p2)
	p12[Pm<kap]=NA
	p12[!is.na(p12)]=1
	image(p1,p2,p12,col="grey90",pch=16,cex=.1,add=T)
	legend("topright",c("Dual DL","Estimated Wald","Null Wald"),lty=1:3,lwd=2,cex=1.5,bg="grey96")
	q=qchisq(lambda,df=2)
	p1.hat=m1.obs/n1;p2.hat=m2.obs/n2
	DIS=p1.hat*(1-p1.hat)/n1*q;DIS[DIS<0]=0
	sqq=sqrt(DIS)
	p1=seq(from=p1.hat-sqq,to=p1.hat+sqq,length=200)
	DIS=p2.hat*(1-p2.hat)*(q-n1*(p1-p1.hat)^2/p1.hat/(1-p1.hat))/n2;DIS[DIS<0]=0
	sqq=sqrt(DIS)
	p2.up=p2.hat+sqq
	p2.low=p2.hat-sqq
	lines(p1,p2.up,lty=2,lwd=2)
	lines(p1,p2.low,lty=2,lwd=2)	
	points(p1.hat,p2.hat,pch=16,cex=1.5)
	
	DIS=q^2/n1^2+4*p1.hat*(1-p1.hat)*q/n1;DIS[DIS<0]=0
	x1=(2*p1.hat+q/n1-sqrt(DIS))/2/(1+q/n1)
	x2=(2*p1.hat+q/n1+sqrt(DIS))/2/(1+q/n1)
	p1s=seq(from=x1,to=x2,length=N)
	q2=1/n2*(q-n1*(p1s-p1.hat)^2/p1s/(1-p1s))
	DIS=q2^2+4*p2.hat*(1-p2.hat)*q2;DIS[DIS<0]=0
	p2s.low=(2*p2.hat+q2-sqrt(DIS))/2/(1+q2)
	p2s.up=(2*p2.hat+q2+sqrt(DIS))/2/(1+q2)
	lines(p1s,p2s.low,lty=3,lwd=2)
	lines(p1s,p2s.up,lty=3,lwd=2)	
	print(Sys.time()-t0)
}
}
