#Onder Model
b0=1 #Type of investment equation
b1=2 #Type of Model (Labor Constrained (1) or Capital Constrained (2))

b2=1 #a0
b3=0.025 #a_Gr
b4=0.1 #a_Resp_K_Gr
b5=1 #a_Resp_Emp
b6=0.95 #a_Resp_Emp_Normal	
b8=0.9 #Consumption propensity
b10=0.5 #Mark_Up_Rt(0)	
b12=1 #Inv_Resp_To_K_Prod_Grwth
b13=0.06 #Inv_Resp_To_K_Util
b14=3 #K/Q
b15=0.8 #Target Capacity Util Rt
b16=0.1 #dep_rt
b18=400 #K0
b19=100 #L
b21=0.1 #MarkUp_Adjstmt_Parameter
b22=1 #Target_Demand/Supply_Rt	1
b24=0 #HH_Autonomous_Resp_to_Wealth	0
b26=0 #Gov_Autonomous_Spending_Rt	0

# time variable
t0=1

I=Ir=u=K=D=Y=Ku=sigma=Sc=Sl=M=Khat=a=ahat=vector()




I[t0]=b18*b16
Ir[t0]=b16
K[t0]=b18
M[t0]=b10
sigma[t0]=1/(1+M[t0])
D[t0]=I[t0]/(1-b8*sigma[t0])
Y[t0]=D[t0]
Ku[t0]=b14*Y[t0]
u[t0]=Ku[t0]/K[t0]
M[2]=M[t0]
Sc[t0]=K[t0]/b14
Sl[t0]=a[t0]*b19  
a[t0]=b2
Khat[t0]=0

Sl
  
for(t in 2:100) {
  K[t]=K[t-1]*(1-b16+Ir[t-1])  
  u[t]=Ku[t-1]/K[t-1]
  
  if (t!=2) {
    if (b1==2) {M[t]=M[t-1]*(1+b21*(D[t-1]/Sc[t-1]-b22))} 
    else {M[t]=M[t-1]*(1+b21*(D[t-1]/Sl[t-1]-b22))}
  }
  
  Khat[t]=(K[t]-K[t-1])/K[t-1]
  a[t]=a[t-1]*(1+b3+b4*Khat[t-1]+b5*(Y[t-1]/(a[t-1]*b19)-b6))
  ahat[t]=(a[t]-a[t-1])/a[t-1]
    Ir[t]=b16+b13*(u[t]-b15)+b12*b3
  I[t]=Ir[t]*K[t]
  Sc[t]=K[t]/b14
  Sl[t]=a[t]*b19  
  sigma[t]=1/(1+M[t])
  D[t]=I[t]/(1-b8*sigma[t])
  Y[t]=D[t]
  Ku[t]=b14*Y[t]

}

par(mfrow=c(1,2), mar=c(3.1, 2.9, 2.5, 1), mgp=c(2, 1, 0), las=0, cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
plot( K, type="l", main="Total Capital", xlab="t", ylab="K",col="blue" )
plot( Ku, type="l", main="Used Capital", xlab="t", ylab="Ku",col="blue" )
plot( u, type="l", main="Utilization Capacity", xlab="t", ylab="u",col="blue" )
plot( M, type="l", main="Mark-up", xlab="t", ylab="M",col="blue" )
plot( sigma, type="l", main="w/a", xlab="t", ylab="w/a",col="blue" )
plot( Ir, type="l", main="Investment Rate", xlab="t", ylab="Ir",col="blue" )
plot( I, type="l", main="Investment", xlab="t", ylab="I",col="blue" )
plot( Sc, type="l", main="A. Supply (Capital)", xlab="t", ylab="Sc",col="blue" )
plot( Sl, type="l", main="A. Supply (Labor)", xlab="t", ylab="Sl",col="blue" )
plot( D, type="l", main="Aggregate Demand", xlab="t", ylab="D",col="blue" )
plot( D/Sc, type="l", main="Demand/C Supply Ratio", xlab="t", ylab="D/S",col="blue" )
plot( D/Sl, type="l", main="Demand/L Supply Ratio", xlab="t", ylab="D/S",col="blue" )
plot( a, type="l", main="Productivity", xlab="t", ylab="a",col="blue" )
plot( ahat, type="l", main="Productivity Growth", xlab="t", ylab="a",col="blue" )









#("Results.pdf") 
# Close the pdf file
#dev.off() 




