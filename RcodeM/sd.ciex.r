sd.ciex <-
function(alpha=0.05)
{
dump("sd.ciex","c:\\Projects\\Mode\\sd.ciex.r")
Y=c(2.3,4.1,3.5,5.2,1.4)
n=length(Y)
S=var(Y)*(n-1)
qL_ET=qchisq(alpha/2,df=n-1)
qU_ET=qchisq(1-alpha/2,df=n-1)
ciL_ET=sqrt(S/qU_ET);ciU_ET=sqrt(S/qL_ET)
qLU_SHORT=var.ql(n=n,adj=0,alpha=alpha)
ciL_SHORT=sqrt(S/qLU_SHORT[2]);ciU_SHORT=sqrt(S/qLU_SHORT[1])
out=as.data.frame(cbind(c(ciL_ET,ciU_ET),cbind(c(ciL_SHORT,ciU_SHORT))))
names(out)=c("Equal-tail CI","Short CI")
row.names(out)=c("Lower limit","Upper limit")
print(out)
L_ET=ciU_ET-ciL_ET
L_SHORT=ciU_SHORT-ciL_SHORT
cat("% CI sigma short =",round((L_ET-L_SHORT)/L_SHORT*100),"\n")
}
