pbinCI <-
function(n=10,m=2,lambda=0.95,maxit=100,eps=10^-6)
{
dump("pbinCI","c:\\Projects\\Mode\\pbinCI.r")
alpha=1-lambda
z1a=qnorm(1-alpha/2)
p.hat=m/n
#Wald
pL.WALD=p.hat-z1a*sqrt(p.hat*(1-p.hat)/n)
pU.WALD=p.hat+z1a*sqrt(p.hat*(1-p.hat)/n)
# Agresti-Coull (equal-tail)
nw=n+z1a^2
pw=(m+z1a^2)/nw
pL.AC=p.hat-z1a*sqrt(pw*(1-pw)/nw)
pU.AC=p.hat+z1a*sqrt(pw*(1-pw)/nw)
# Clopper-Pearson				
pL.CP=qbeta(alpha/2,shape1=m+.5,shape2=n-m+1)
pU.CP=qbeta(1-alpha/2,shape1=m+1,shape2=n-m)
#Beta equal-tail
pL.ETBETA=qbeta((1-lambda)/2,m+1,n-m+1)
pU.ETBETA=qbeta((1+lambda)/2,m+1,n-m+1)
#Optimal MCL
pL.MCL=pL.ETBETA;pU.MCL=pU.ETBETA
for(it in 1:maxit)
{
	LHS1=pbeta(pU.MCL,m+1,n-m+1)-pbeta(pL.MCL,m+1,n-m+1)-lambda
	LHS2=(m+1)*log(pU.MCL/pL.MCL)-(n-m+1)*log((1-pL.MCL)/(1-pU.MCL))
	H11=-dbeta(pL.MCL,m+1,n-m+1);H12=dbeta(pU.MCL,m+1,n-m+1)
	H21=-(m+1)/pL.MCL+(n-m+1)/(1-pL.MCL);H22=(m+1)/pU.MCL-(n-m+1)/(1-pU.MCL)
	detH=H11*H22-H12*H21
	del1=LHS1*H22-LHS2*H12
	del2=LHS2*H11-LHS1*H21
	ma=max(abs(c(LHS1,LHS2)))
	#print(c(it,pL.MCL,pU.MCL,ma))
	if(ma<eps) break
	pL.MCL=pL.MCL-del1/detH
	pU.MCL=pU.MCL-del2/detH					
}
out=as.data.frame(cbind(c(pL.WALD,pU.WALD),c(pL.AC,pU.AC),c(pL.CP,pU.CP),c(pL.ETBETA,pU.ETBETA),c(pL.MCL,pU.MCL)))
names(out)=c("Asymptotic/Wald","Agresti-Coul","Clopper-Pearson","MET","MCL2")
row.names(out)=c("Lower","Upper")
cat("Estimation of binomial probability n=",n,", m=",m,"\n",lambda*100,"% Confidence intervals\n",sep="")
out
}
