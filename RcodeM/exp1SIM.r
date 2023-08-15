exp1SIM.out <-
structure(c(0.10000000000000001, 0.20000000000000001, 0.30000000000000004, 
0.40000000000000002, 0.5, 0.59999999999999998, 0.70000000000000007, 
0.80000000000000004, 0.90000000000000002, 1, 1.1000000000000001, 
1.2000000000000002, 1.3000000000000003, 1.4000000000000001, 1.5, 
0.89148000000000005, 0.89140457022851138, 0.88687981918372649, 
0.88435836901284459, 0.88226464937560034, 0.8771300987204389, 
0.87361470169742883, 0.87067271049781381, 0.86571029885657491, 
0.86094772720417456, 0.85618416154380173, 0.85251236788279983, 
0.84791752826896116, 0.84251650548163703, 0.84348265749823037, 
0.94940000000000002, 0.94886744337216866, 0.94482503425308273, 
0.94108879196510742, 0.9383005123278898, 0.93512084743386936, 
0.93228321208841858, 0.92900838381002049, 0.92572254625400296, 
0.92092218904272116, 0.91776921451640492, 0.91414523068242504, 
0.91129447128720864, 0.90548971309737725, 0.90547072504803316, 
0.94954598416971348, 0.95131014481077925, 0.94872363145487881, 
0.94910285213431522, 0.94961252028096055, 0.95051352932867894, 
0.94977196593615043, 0.95056404455800514, 0.95095354078446037, 
0.95100053645353599, 0.95064142836876553, 0.95064808616807639, 
0.95071795600869269, 0.95008535891390944, 0.95082401083072909
), .Dim = c(15L, 4L), .Dimnames = list(NULL, c("s2s", "cpZ", 
"cpT", "cpEX")))
exp1SIM <-
function(job=1,alpha=0.05,nSim=100000,ss=3)
{
dump(c("exp1SIM.out","exp1SIM"),"c:\\Projects\\Mode\\exp1SIM.r")
set.seed(ss)
if(job==1)
{
	t0=Sys.time()
	#install.packages("nleqslv")
	library(nleqslv)
	Qbf=function(b,xd,y,tnm)
	{
		n=length(xd)
		f=exp(b*xd)
		g=xd*f
		res=y-f
		sden=max(sum(res^2)*sum(g^2)-(sum(g*res))^2,10^-10)
		Qb=sum(g*res)*sqrt(n-1)/sqrt(sden)-tnm
		return(Qb)	
	}
	
	s2s=seq(from=.1,to=1.5,by=.1)
	ns2s=length(s2s)
	b.true=.1;x=c(-4,-1,2,3,4,5);n=length(x)	
	tcr=qt(1-alpha/2,df=n-1)
	Zcr=qnorm(1-alpha/2)
	cpZ=cpT=cpEX=rep(NA,ns2s)	
	for(j in 1:ns2s)
	{
		up.NLS=low.NLS=up.T=low.T=up.ex=low.ex=rep(NA,nSim)
		for(i in 1:nSim)
		{
			y=exp(b.true*x)+rnorm(n)*sqrt(s2s[j])
			o=try(nls(y~exp(b*x),start=c(b=b.true)))
			if(attr(o,"class")!="try-error")
			{
				so=summary(o)$coefficients
				seb=so[1,2]
				b.NLS=so[1,1]
				up.NLS[i]=b.NLS+Zcr*seb;low.NLS[i]=b.NLS-Zcr*seb			
				up.T[i]=b.NLS+tcr*seb;low.T[i]=b.NLS-tcr*seb
				for(k in 1:100)
				{
					xtt=b.NLS-0.1*k*seb
					if(Qbf(xtt,xd=x,y=y,tcr)>0) break
				}
			
				out=nleqslv(x=xtt,fn=Qbf,xd=x,y=y,tnm=tcr)
				if(abs(out$fvec)<10^-7)	low.ex[i]=out$x
			
				for(k in 1:100)
				{
					xtt=b.NLS+0.1*k*seb
					if(Qbf(xtt,xd=x,y=y,-tcr)<0) break
				}
				out=nleqslv(x=xtt,fn=Qbf,xd=x,y=y,tnm=-tcr)
				if(abs(out$fvec)<10^-7)	up.ex[i]=out$x
			}
		}
		cpZ[j]=mean(low.NLS<b.true & b.true<up.NLS,na.rm=T)
		cpT[j]=mean(low.T<b.true & b.true<up.T,na.rm=T)
		cpEX[j]=mean(low.ex<b.true & b.true<up.ex,na.rm=T)
	}
	print(Sys.time()-t0)
	#exp1SIM.out=exp1SIM(job=1)
	return(cbind(s2s,cpZ,cpT,cpEX))	
}
if(job==2)
{
	par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=1.4)
	n=length(exp1SIM.out[,1])
	matplot(sqrt(exp1SIM.out[,1]),exp1SIM.out[,c(2,3,4)],ylim=c(.82,1),lwd=2,col=1,lty=1,type="o",pch=1:3,cex=1.5,xlab="Variance",ylab="Coverage probability")
	segments(-1,1-alpha,10,1-alpha)
	text(.5,.98,"n = 6",font=3,cex=2)
	legend("bottomleft",c("N-Wald/asymptotic CI","T-Wald/asymptotic CI","Exact pivotal CI"),pch=1:3,lty=1,lwd=2,cex=1.6,bg="gray95")
	
}

}
