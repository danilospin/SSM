#Onder Model 2020

#Setting the libraries
library(deSolve)
library(rootSolve)
library(psych)
library(scatterplot3d)

#SETTING THE MODEL
#print results in a .pdf file? 1 - Yes, 2 - No
print_results=0

#Type of investment equation (Using Labor Pdtv growth (1) or fixed Labor Pdtv (2))
b0=2 
#Type of Model (Labor Constrained (1) or Capital Constrained (2))
b1=2 

#PARAMETERS
b6=0.05 #Inv Response to Unemployment	
b7=0.95 #Natural full employment
b8=0.9 #Consumption propensity
b10=0.5 #Mark_Up_Rt(0)
b12=0.05 #Inv_Resp_To_K_Prod_Grwth
b13=0.06 #Inv_Resp_To_K_Util
b14=3 #K/Q
b15=1 #Target Capacity Util Rt
b16=0.1 #dep_rt
b18=300 #K0
b19=100 #L
b21=0.1 #MarkUp_Adjstmt_Parameter
b22=1 #Target_Demand/Supply_Rt	1
b24=0 #HH_Autonomous_Resp_to_Wealth	0
b26=0 #Gov_Autonomous_Spending_Rt	0

a1=1
a2=2
a3=4
a4=8
# Declareing vectors
t0=1
I=Ir=u=K=D=Y=Ku=sigma=Sc=Sl=M=Khat=a=ahat=e=K1=K2=K3=K4=Ir1=Ir2=Ir3=Ir4=I1=I2=I3=I4=abar=vector()

#Initial conditions
K1[t0]=b18
K2[t0]=0
K3[t0]=0
K4[t0]=0
Ir1[t0]=0.1
Ir2[t0]=0.1
Ir3[t0]=0.1
Ir4[t0]=0.1
I1[t0]=Ir1[t0]*K1[t0]
I2[t0]=Ir2[t0]*K2[t0]
I3[t0]=Ir3[t0]*K3[t0]
I4[t0]=Ir4[t0]*K4[t0]

I[t0]=b18*b16 #Investment
Ir[t0]=b16 #Investment Rate
K[t0]=b18 #Actual Capital Stock

M[t0]=b10 #Mark-up
sigma[t0]=1/(1+M[t0]) #Wage Share in total income
D[t0]=I[t0]/(1-b8*sigma[t0]) #Aggregate Demand
Y[t0]=D[t0] #Total Output
Ku[t0]=b14*Y[t0] #Capital Stock at full capacity
u[t0]=Ku[t0]/K[t0] #Utilization capacity
M[2]=M[t0] 
Sc[t0]=K[t0]/b14 #Total Supply with capital constraints  
Sl[t0]=a[t0]*b19 #Total Supply with labor constraints
Khat[t0]=0 #Capital stock growth
abar[t0]=1
ahat[t0]=0 #Labor productivity growth
e[t0]=Y[t0]/(b19*abar[t0]) #Total Employment


# Dynamic Equations  
for(t in 2:1000) {
  
  u[t]=Ku[t-1]/K[t-1]
  
  K1[t]=K1[t-1]*(1-b16+Ir1[t-1])
  K2[t]=K2[t-1]*(1-b16+Ir2[t-1])
  K3[t]=K3[t-1]*(1-b16+Ir3[t-1])
  K4[t]=K4[t-1]*(1-b16+Ir4[t-1])
  
  #Shock
  if(t==128)
  { K1[t]=K1[t-1]*(1-b16+Ir1[t-1])*99/100
    K2[t]=K1[t-1]*(1-b16+Ir1[t-1])/100}
  if(t==403)
  { K2[t]=K2[t-1]*(1-b16+Ir2[t-1])*99/100
  K3[t]=K2[t-1]*(1-b16+Ir2[t-1])/100}
  if(t==603)
  { K3[t]=K3[t-1]*(1-b16+Ir3[t-1])*99/100
  K4[t]=K3[t-1]*(1-b16+Ir3[t-1])/100}

  K[t]=K1[t]+K2[t]+K3[t]+K4[t]  
  
  Ir1[t]=b16+b13*(u[t]-b15)+b12*(a1-abar[t-1])+b6*(b7-e[t-1])
  Ir2[t]=b16+b13*(u[t]-b15)+b12*(a2-abar[t-1])+b6*(b7-e[t-1])
  Ir3[t]=b16+b13*(u[t]-b15)+b12*(a3-abar[t-1])+b6*(b7-e[t-1])
  Ir4[t]=b16+b13*(u[t]-b15)+b12*(a4-abar[t-1])+b6*(b7-e[t-1])
  
  I1[t]=Ir1[t]*K1[t]
  I2[t]=Ir2[t]*K2[t]
  I3[t]=Ir3[t]*K3[t]
  I4[t]=Ir4[t]*K4[t]
  
  I[t]=I1[t]+I2[t]+I3[t]+I4[t]
  abar[t]=(a1*K1[t]+a2*K2[t]+a3*K3[t]+a4*K4[t])/K[t]
  
  M[t]=M[t-1]*(1+b21*(D[t-1]/Sc[t-1]-b22))
  M[2]=M[t0]
  
  Khat[t]=(K[t]-K[t-1])/K[t-1]
  ahat[t]=(abar[t]-abar[t-1])/abar[t-1]
  Sc[t]=K[t]/b14
  Sl[t]=abar[t]*b19  
  sigma[t]=1/(1+M[t])
  D[t]=I[t]/(1-b8*sigma[t])
  Y[t]=D[t]
  Ku[t]=b14*Y[t]
  e[t]=Y[t]/(b19*abar[t])
}

#Graphs
par(mfrow=c(1,2), mar=c(3.1, 2.9, 2.5, 1), mgp=c(2, 1, 0), las=0, cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
plot( M, type="l", main="Mark-up", xlab="t", ylab="M",col="blue" )
plot( sigma, type="l", main="w/a", xlab="t", ylab="w/a",col="blue" )
plot( abar, type="l", main="Productivity", xlab="t", ylab="a",col="blue" )
plot( ahat, type="l", main="Productivity Growth", xlab="t", ylab="ahat",col="blue" )

plot( K1, type="l", main="Capital 1", xlab="t", ylab="K",col="blue" )
plot( K2, type="l", main="Capital 2", xlab="t", ylab="K",col="blue" )
plot( K3, type="l", main="Capital 3", xlab="t", ylab="K",col="blue" )
plot( K4, type="l", main="Capital 4", xlab="t", ylab="K",col="blue" )
plot( K, type="l", main="Capital Total", xlab="t", ylab="K",col="blue" )
plot( Khat, type="l", main="Capital Growth Rate", xlab="t", ylab="Khat",col="blue" )
plot( Ku, type="l", main="Used Capital", xlab="t", ylab="Ku",col="blue" )
plot( u, type="l", main="Utilization Employment", xlab="t", ylab="u(blue) e(red)",col="blue" )
lines( e, type="l", xlab="t", ylab="e",col="red" )

plot( I1, type="l", main="Investment 1", xlab="t", ylab="I",col="blue" )
plot( I2, type="l", main="Investment 2", xlab="t", ylab="I",col="blue" )
plot( I3, type="l", main="Investment 3", xlab="t", ylab="I",col="blue" )
plot( I4, type="l", main="Investment 4", xlab="t", ylab="I",col="blue" )

plot( Sc, type="l", main="A. Supply (Capital)", xlab="t", ylab="Sc",col="blue" )
plot( D, type="l", main="Aggregate Demand", xlab="t", ylab="D",col="blue" )
plot( D/Sc, type="l", main="Demand/C Supply Ratio", xlab="t", ylab="D/S",col="blue" )
plot( D/Sl, type="l", main="Demand/L Supply Ratio", xlab="t", ylab="D/S",col="blue" )



# Printing the results in a .pdf file
if (print_results==1)
{
pdf("Results_Onder_Model_vintage.pdf") 
  # Graph Characteristics
  # par(mfrow=c(1,1), mar=c(3.1, 2.9, 2.5, 1), mgp=c(2, 1, 0), las=0, cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), lwd = 2, mgp=c(3, 1, 0), las=0, cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
  
  #Graphs
  plot( M, type="l", main="Mark-up", xlab="t", ylab="M",col="blue" )
  plot( sigma, type="l", main="w/a", xlab="t", ylab="w/a",col="blue" )
  plot( abar, type="l", main="Productivity", xlab="t", ylab="a",col="blue" )
  plot( ahat, type="l", main="Productivity Growth", xlab="t", ylab="ahat",col="blue" )
  
  plot( K1, type="l", main="Capital 1", xlab="t", ylab="K",col="blue" )
  plot( K2, type="l", main="Capital 2", xlab="t", ylab="K",col="blue" )
  plot( K3, type="l", main="Capital 3", xlab="t", ylab="K",col="blue" )
  plot( K4, type="l", main="Capital 4", xlab="t", ylab="K",col="blue" )
  plot( K, type="l", main="Capital Total", xlab="t", ylab="K",col="blue" )
  plot( Khat, type="l", main="Capital Growth Rate", xlab="t", ylab="Khat",col="blue" )
  plot( Ku, type="l", main="Used Capital", xlab="t", ylab="Ku",col="blue" )
  plot( u, type="l", main="Utilization Employment", xlab="t", ylab="u(blue) e(red)",col="blue" )
  lines( e, type="l", xlab="t", ylab="e",col="red" )
  
  plot( I1, type="l", main="Investment 1", xlab="t", ylab="I",col="blue" )
  plot( I2, type="l", main="Investment 2", xlab="t", ylab="I",col="blue" )
  plot( I3, type="l", main="Investment 3", xlab="t", ylab="I",col="blue" )
  plot( I4, type="l", main="Investment 4", xlab="t", ylab="I",col="blue" )
  
  plot( Sc, type="l", main="A. Supply (Capital)", xlab="t", ylab="Sc",col="blue" )
  plot( D, type="l", main="Aggregate Demand", xlab="t", ylab="D",col="blue" )
  plot( D/Sc, type="l", main="Demand/C Supply Ratio", xlab="t", ylab="D/S",col="blue" )
  plot( D/Sl, type="l", main="Demand/L Supply Ratio", xlab="t", ylab="D/S",col="blue" )
#Close the pdf file
dev.off() 
}

