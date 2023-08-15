exp2 <-
function(job=0,n1=4,l1=.3,n2=10,l2=2,alpha=0.05,N=100,nSim=1000000,ss=2)
{
dump("exp2","c:\\Projects\\Mode\\exp2.r")
q12=function(n1,n2,d=0,alpha=0.05,maxit=1000,eps=0.001)
{
	al0=alpha/2
	q1=qf(al0,df1=2*n1,df2=2*n2)
	q2=qf(1-al0,df1=2*n1,df2=2*n2)
	H=matrix(ncol=2,nrow=2);b=rep(NA,2)
	for(it in 1:maxit)
	{
		b[1]=pf(q2,df1=2*n1,df2=2*n2)-pf(q1,df1=2*n1,df2=2*n2)-(1-alpha)
		b[2]=(n1+d)*log(q2/q1)-(n1+n2)*log((n2+n1*q2)/(n2+n1*q1))
		#print(c(it,q1,q2,b))
		H[1,1]=-df(q1,df1=2*n1,df2=2*n2)
		H[1,2]=df(q2,df1=2*n1,df2=2*n2)
		H[2,1]=-(n1+d)/q1+n1*(n1+n2)/(n2+n1*q1)
		H[2,2]=(n1+d)/q2-n1*(n1+n2)/(n2+n1*q2)
		delta=.25*solve(H)%*%b
		if(max(abs(delta))<eps) break
		q1=q1-delta[1];q2=q2-delta[2]	
	}
	return(c(q1,q2))
}
if(job==0) # cheking X_bar/(nu*Y_bar) ~ F(2*n1,2*n2)
{
	l1=runif(1);l2=runif(1)
	Xi=matrix(rexp(n1*nSim,rate=l1),ncol=n1)
	S1=rowMeans(Xi)
	Yj=matrix(rexp(n2*nSim,rate=l2),ncol=n2)
	S2=rowMeans(Yj)
	nu=l2/l1
	s=S1/S2
	s=sort(s)
	plot(s,(1:nSim)/nSim,type="s")
	x=seq(from=0,to=2*n1,length=1000)
	lines(x,pf(x/nu,df1=2*n1,df2=2*n2),col=2)
}
if(job==1) #Power functions
{
	n1=1;l1=.3;n2=10;l2=.3
	nu0=l2/l1
	nu=seq(from=.1*nu0,to=3*nu0,length=N)
	q1=qf(alpha/2,df1=2*n1,df2=2*n2)
	q2=qf(1-alpha/2,df1=2*n1,df2=2*n2)
	Pnu0.eq=1-pf(q2*nu0/nu,df1=2*n1,df2=2*n2)+pf(q1*nu0/nu,df1=2*n1,df2=2*n2)
	q12.unb=q12(n1=n1,n2=n2)
	Pnu0.unb=1-pf(q12.unb[2]*nu0/nu,df1=2*n1,df2=2*n2)+pf(q12.unb[1]*nu0/nu,df1=2*n1,df2=2*n2)
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	matplot(nu,cbind(Pnu0.eq,Pnu0.unb),lwd=2,ylim=c(0,.25),col=1,type="l",xlab="",ylab="Power")
	segments(-1,alpha,1000,alpha,lty=2)
		
	#check through simulations
	nuSIM=nu[seq(from=1,to=N,length=10)];N5=length(nuSIM)
	powSIM.eq=powSIM.unb=rep(NA,N5)
	for(i in 1:N5)
	{
		Xi=matrix(rexp(n1*nSim,rate=1),ncol=n1)
		S1=rowMeans(Xi)
		Yj=matrix(rexp(n2*nSim,rate=nuSIM[i]),ncol=n2)
		S2=rowMeans(Yj)
		powSIM.eq[i]=1-mean(S1/S2/nu0<q2 & S1/S2/nu0>q1)		
	}	
	
	for(i in 1:N5)
	{
		Xi=matrix(rexp(n1*nSim,rate=1),ncol=n1)
		S1=rowMeans(Xi)
		Yj=matrix(rexp(n2*nSim,rate=nuSIM[i]),ncol=n2)
		S2=rowMeans(Yj)
		powSIM.unb[i]=1-mean(S1/S2/nu0<q12.unb[2] & S1/S2/nu0>q12.unb[1])			
	}
	points(nuSIM,powSIM.eq,cex=1.5)
	points(nuSIM,powSIM.unb,cex=1.5,pch=2)	
	segments(nu0,-1,nu0,2,lty=2)
	legend("topleft",c("Equal-tail","Unbiased"),lty=1:2,pch=1:2,lwd=2,cex=1.8,bg="grey96")
	
}
if(job==2) #optimal design
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	n1=4;l1=.3;n2=10;l2=.3
	nu0=l1/l1;nu.alt=2
	n1=seq(from=20,to=80,length=N);n2=seq(from=20,to=80,length=N)
	pow=matrix(nrow=N,ncol=N)
	for(i1 in 1:N)
	for(i2 in 1:N)
	{
		q12.unb=q12(n1=n1[i1],n2=n2[i2])
		pow[i1,i2]=pf(q12.unb[2]*nu0/nu.alt,df1=2*n1[i1],df2=2*n2[i2],lower.tail=F)+pf(q12.unb[1]*nu0/nu.alt,df1=2*n1[i1],df2=2*n2[i2])
	}
	contour(n1,n2,pow,levels=c(.8,.9),lwd=2,xlab="Sample size, n1",ylab="Sample size, n2",labcex=1)
	text(60,65,paste("nu0 =",nu0,"\nnu.alt =",nu.alt),cex=1.5,font=3)
	for(lev in c(.8,.9))
	{
		con=contourLines(n1,n2,levels=lev,pow)
		n1c=con[[1]]$x;n2c=con[[1]]$y
		N=length(n1)
		n11=n1c[2:N];n22=n2c[2:N]
		sl=abs((n22-n2c[1:(N-1)])/(n11-n1c[1:(N-1)])+1)
		n1.opt=n1c[2:N][sl==min(sl)]
		n2.opt=n2c[2:N][sl==min(sl)]
		points(n1.opt,n2.opt,pch=1,cex=2)
		text(n1.opt+3,n2.opt+.5,paste("n1 = ",round(n1.opt),", n2 = ",round(n2.opt),sep=""),adj=0,cex=1.25)
		segments(n1.opt-10,n2.opt+10,n1.opt+10,n2.opt-10,lty=2)
	}	
}
if(job==3) #short CI
{
	t0=Sys.time()
	n1=2;n2=10
	nuSIM=seq(from=.2,to=2,length=N)
	q1=qf(alpha/2,df1=2*n1,df2=2*n2)
	q2=qf(1-alpha/2,df1=2*n1,df2=2*n2)
	covpr.et=covpr.short=wid.et=wid.sh=rep(NA,N)
	for(i in 1:N)
	{
		Xi=matrix(rexp(n1*nSim,rate=1),ncol=n1)
		S1=rowMeans(Xi)
		Yj=matrix(rexp(n2*nSim,rate=nuSIM[i]),ncol=n2)
		S2=rowMeans(Yj)
		S=S1/S2
		L.et=S/q2;U.et=S/q1
		covpr.et[i]=mean(nuSIM[i]>L.et & nuSIM[i]<U.et)
		wid.et[i]=mean((U.et-L.et)/nuSIM[i])
		LU.q=q12(n1=n1,n2=n2,d=1)
		L.short=S/LU.q[2];U.short=S/LU.q[1]
		covpr.short[i]=mean(nuSIM[i]>L.short & nuSIM[i]<U.short)				
		wid.sh[i]=mean((U.short-L.short)/nuSIM[i])
	}
	par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),cex.lab=1.5)
	matplot(nuSIM,cbind(covpr.et,covpr.short),ylim=c(0.9,1),col=1,lty=1:2,lwd=2,xlab="True nu",ylab="Probability",main="Coverage probability",type="l")
	text(1,.92,paste("n1 = ",n1,", n2 = ",n2,sep=""),cex=1.5)
	legend("topleft",c("Equal-tail CI","Short CI"),col=1,lty=1:2,lwd=2,cex=1.5,bg="grey96")
	matplot(nuSIM,cbind(wid.et,wid.sh),ylim=c(5,10),col=1,lty=1:2,lwd=2,xlab="True nu",ylab="Length",main="Relative length of CI",type="l")
	print(Sys.time()-t0)
}
if(job==4) #CI for nu on the log scale
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	nn=matrix(c(2,10,4,5),nrow=2)
	nuSIM=seq(from=.2,to=2,length=N)
	wid.et=wid.sh=rep(NA,N)
	plot(nuSIM,nuSIM,ylim=c(2.5,5),type="n")
	for(inn in 1:2)
	{
		n1=nn[1,inn];n2=nn[2,inn]
		q1=qf(alpha/2,df1=2*n1,df2=2*n2)
		q2=qf(1-alpha/2,df1=2*n1,df2=2*n2)
		for(i in 1:N)
		{
			Xi=matrix(rexp(n1*nSim,rate=1),ncol=n1)
			S1=rowMeans(Xi)
			Yj=matrix(rexp(n2*nSim,rate=nuSIM[i]),ncol=n2)
			S2=rowMeans(Yj)
			S=S1/S2
			L.et=S/q2;U.et=S/q1
			wid.et[i]=mean(log(U.et/L.et))
			LU.q=q12(n1=n1,n2=n2,d=2)
			L.short=S/LU.q[2];U.short=S/LU.q[1]
			wid.sh[i]=mean(log(U.short/L.short))
		}
		lines(nuSIM,wid.et,lty=1,lwd=2);lines(nuSIM,wid.sh,lty=2,lwd=2)
		print(cbind(nuSIM,wid.et,wid.sh))		
	}	
}
if(job==4.1) #Demonstartion of MC estimator and shortest CI
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.5)
	n1=2;n2=10
	S=1
	al=ciL=ciU=seq(from=0.05,to=0.99,length=200);nal=length(al)
	adj=1
	for(i in 1:nal)
	{		
		qLU=q12(n1=n1,n2=n2,d=adj,alpha=al[i])
		ciL[nal-i+1]=S/qLU[2];ciU[nal-i+1]=S/qLU[1]
	}
	matplot(al,cbind(ciL,ciU),lwd=2,xlim=c(0.05,1),type="l",lty=1,col=1,xlab="Confidence level",ylab="Lower and upper limits")	
	lim=S*n1*(n2-1)/n2/(n1+1)
	segments(-1,lim,2,lim,lty=2,lwd=2)	
	legend("topleft",c("Confidence limits","MC estimator of nu"),lty=1:2,lwd=2,cex=1.5,bg="grey95")
}

}
