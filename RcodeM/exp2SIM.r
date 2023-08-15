exp2SIM <-
function(job=1,alpha=0.05,nSim=10000,ss=3)
{
dump(c("exp2SIM","exp2SIM_out"),"c:\\Projects\\Mode\\exp2SIM.r")
set.seed(ss)
t0=Sys.time()
#install.packages("numDeriv")
library(numDeriv)
cieq=function(b,xdat,y,q1a) #Profile equation for beta
{
	f0=exp(b*xdat);g0=f0*xdat
	A=sum(y*f0);B=sum(y*g0);C=sum(f0*g0);D=sum(f0^2)
	gh=sum(y*f0)/D
	gb=(A*g0+B*f0)/D-2*A*C/D^2*f0
	r=y-gh*f0
	num=sum(gb*r)*sqrt(n-2)
	sqden=sum(r^2)*sum(gb^2)-sum(gb*r)^2
	return(num/sqrt(sqden)-q1a)	
}
	
ciEx=function(b,xdat,y,q1a) #Solving the profile equation for beta
{
	for(iter in 1:100)
	{
		Q0=try(cieq(b=b,xdat=xdat,y=y,q1a=q1a),silent=T)
		if(!is.numeric(Q0)) return(NA)
		der=try(grad(func=cieq,x=b,xdat=xdat,y=y,q1a=q1a),silent=T)
		if(!is.numeric(der)) return(NA)			
		b.new=b-Q0/der
		if(abs(Q0)<0.0001) return(b)
		b=b.new
		#print(c(iter,b,Q0))
	}
	return(NA)	
}
#doubly noncentral t-density
prtf=function(v,nm,q,d1,d2) dchisq(v,df=nm,ncp=d2)*pnorm(q/sqrt(nm)*sqrt(v)-d1)

if(job==1) #simulations and saving them in file exp2SIM_out
{	
	b.true=.2
	x=c(-4,-1,2,3,4,5)
	g.true=1
	n=length(x)
	q1a=qt(1-alpha/2,n-2)
	Z1a=qnorm(1-alpha/2)
	N=20
	sigs=seq(from=0.1,to=.5,length=N)
	cpNLS=cpEP=rep(NA,N)
	#CIs
	for(iss in 1:N)
	{
		bs=lhs=seq(from=0,to=2,length=N)
		low=up=lowNLS=upNLS=rep(NA,nSim)
		for(i in 1:nSim)
		{
			y=g.true*exp(b.true*x)+rnorm(n)*sigs[iss]
			o=try(nls(y~g*exp(b*x),start=c(g=g.true,b=b.true)),silent=T)
			if(attr(o,"class")!="try-error")
			{
				so=summary(o)$coefficients
				bnls=so[2,1]
				se=so[2,2]
				lowNLS[i]=bnls-Z1a*se;upNLS[i]=bnls+Z1a*se
				low[i]=ciEx(b=bnls-se*2,xdat=x,y=y,q1a=q1a)
				if(is.na(low[i])) low[i]=-1000
				up[i]=ciEx(b=bnls+se*2,xdat=x,y=y,q1a=-q1a)
				if(is.na(up[i])) up[i]=1000			
			}			
		}
		cpNLS[iss]=mean(lowNLS<b.true & upNLS>b.true,na.rm=T)
		cpEP[iss]=mean(low<b.true & up>b.true,na.rm=T)
	}
	#Powers
	NBS=20
	bs=seq(from=.3,to=.6,length=NBS)
	powNLS=powEX=rep(NA,NBS)
	b.true=0.5
	sigsPOW=.4
	for(j in 1:NBS)
	{
		Tex=Z=rep(NA,nSim)
		for(i in 1:nSim)
		{
			y=g.true*exp(bs[j]*x)+rnorm(n)*sigsPOW
			o=try(nls(y~g*exp(b*x),start=c(g=g.true,b=bs[j])),silent=T)
			if(attr(o,"class")!="try-error")
			{
				so=summary(o)$coefficients
				bnls=so[2,1]
				se=so[2,2]
				Z[i]=(bnls-b.true)/se
				Tex[i]=cieq(b=b.true,xdat=x,y=y,q1a=0)				
			}
			powNLS[j]=mean(abs(Z)>Z1a,na.rm=T)
			powEX[j]=mean(abs(Tex)>q1a,na.rm=T)
		}	
	}	
	print(Sys.time()-t0)
	#exp2SIM_out=exp2SIM(job=1)
	return(list(cbind(sigs,cpNLS,cpEP),cbind(bs,powNLS,powEX)))
}
if(job==2)#plotting the figure
{
	par(mfrow=c(1,2),mar=c(4.5,4.5,3,1),cex.lab=1.5,cex.main=1.5)
	cp=(exp2SIM_out[[1]])[,2:3]
	N=20
	sigs=seq(from=0.1,to=.5,length=N)
	matplot(sigs,cp,ylim=c(.8,1),col=1,lty=1,pch=1:2,cex=1.5,type="o",ylab="Coverage probability",xlab="sigma",main="Confidence interval")
	segments(-1,1-alpha,10,1-alpha)
	legend("bottomleft",c("Asymptotic normal CI","Pivotal CI"),pch=c(1,2),lty=1,cex=1.5,bg="gray95")
	
	bs=(exp2SIM_out[[2]])[,1]
	cp=(exp2SIM_out[[2]])[,2:3]
	matplot(bs,cp,ylim=c(0,1),col=1,lty=1,pch=1:2,cex=1.5,type="o",ylab="Power",xlab="Alternative beta",main="Power function")
	legend("topleft",c("Z-test","Profile pivotal test","Analytic power"),pch=c(1,2,NA),lty=1,lwd=c(1,1,2),cex=1.6,bg="gray95")	
	b.true=0.5
	segments(b.true,-1,b.true,.6,lty=2)
	segments(-1,alpha,100,alpha,lty=3)
	bss=pow.an=seq(from=0.3,to=0.6,length=100)
	g.true=1
	n=6
	sigsPOW=.4
	x=c(-4,-1,2,3,4,5)
	t1a=qt(1-alpha/2,df=n-2)	
	for(i in 1:100)
	{
		A=sum(exp((b.true+bss[i])*x))
		B=sum(exp(2*b.true*x))
		C=exp(b.true*x)
		D=sum(x*exp(2*b.true*x))
		g0=g.true*A/B
		dff=g.true*exp(bss[i]*x)-g0*exp(b.true*x)
		g=g0*x*exp(b.true*x)-(g0*D-sum(x*exp(b.true*x)*dff))*exp(b.true*x)/sum(C^2)		
		delta1=sum(dff*g)/sqrt(sum(g^2))/sigsPOW
		delta2=sum(dff^2)/sigsPOW^2-sum(dff*g)^2/sum(g^2)/sigsPOW^2				
		pow.an[i]=1-integrate(prtf,nm=n-2,q=t1a,d1=delta1,d2=delta2,lower=0,upper=Inf)$value+integrate(prtf,nm=n-2,q=-t1a,d1=delta1,d2=delta2,lower=0,upper=Inf)$value				
		
	}
	lines(bss,pow.an,lwd=2)	
}
	
}
exp2SIM_out <-
list(structure(c(0.10000000000000001, 0.12105263157894737, 0.14210526315789473, 
0.16315789473684211, 0.18421052631578949, 0.20526315789473684, 
0.22631578947368422, 0.24736842105263157, 0.26842105263157895, 
0.28947368421052633, 0.31052631578947365, 0.33157894736842108, 
0.35263157894736841, 0.37368421052631584, 0.39473684210526316, 
0.41578947368421049, 0.43684210526315792, 0.45789473684210524, 
0.47894736842105268, 0.5, 0.87490000000000001, 0.87409999999999999, 
0.87749999999999995, 0.87860000000000005, 0.87729999999999997, 
0.87319999999999998, 0.87660000000000005, 0.87939999999999996, 
0.87670000000000003, 0.874, 0.87760000000000005, 0.87090000000000001, 
0.87729999999999997, 0.87458745874587462, 0.87390000000000001, 
0.87260000000000004, 0.87590000000000001, 0.87360000000000004, 
0.86625987796338899, 0.86734693877551017, 0.94799999999999995, 
0.94879999999999998, 0.95230000000000004, 0.9506, 0.95230000000000004, 
0.95169999999999999, 0.95150000000000001, 0.95440000000000003, 
0.95150000000000001, 0.95550000000000002, 0.95320000000000005, 
0.95509999999999995, 0.95630000000000004, 0.95569556955695567, 
0.9546, 0.9577, 0.95940000000000003, 0.95589999999999997, 0.9569870961288387, 
0.95558223289315725), .Dim = c(20L, 3L), .Dimnames = list(NULL, 
    c("sigs", "cpNLS", "cpEP"))), structure(c(0.29999999999999999, 
0.31578947368421051, 0.33157894736842103, 0.34736842105263155, 
0.36315789473684212, 0.37894736842105259, 0.39473684210526316, 
0.41052631578947368, 0.4263157894736842, 0.44210526315789472, 
0.45789473684210524, 0.47368421052631582, 0.48947368421052628, 
0.50526315789473686, 0.52105263157894743, 0.5368421052631579, 
0.55263157894736836, 0.56842105263157894, 0.58421052631578951, 
0.59999999999999998, 0.90510000000000002, 0.88859999999999995, 
0.85140000000000005, 0.8327, 0.78779999999999994, 0.74960000000000004, 
0.6734, 0.60129999999999995, 0.50609999999999999, 0.40910000000000002, 
0.29899999999999999, 0.20880000000000001, 0.14760000000000001, 
0.11990000000000001, 0.1812, 0.32419999999999999, 0.52780000000000005, 
0.75509999999999999, 0.91500000000000004, 0.98180000000000001, 
0.54249999999999998, 0.54510000000000003, 0.52310000000000001, 
0.50870000000000004, 0.47049999999999997, 0.43369999999999997, 
0.37759999999999999, 0.32669999999999999, 0.25950000000000001, 
0.19339999999999999, 0.1389, 0.0877, 0.057599999999999998, 0.051700000000000003, 
0.087400000000000005, 0.1739, 0.33129999999999998, 0.55110000000000003, 
0.75760000000000005, 0.90849999999999997), .Dim = c(20L, 3L), .Dimnames = list(
    NULL, c("bs", "powNLS", "powEX"))))
