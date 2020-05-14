#Onder Model

#Setting the libraries
library(deSolve)
library(rootSolve)
library(psych)
library(scatterplot3d)

#Setting the Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwg("C:/LSD/LSD-7.2-stable-2/Example/SantAnna/K+S-original/R/data")
getwd()

#SETTING THE MODEL
#print results in a .pdf file? 1 - Yes, 2 - No
print_results=1

#Type of investment equation (Using Labor Pdtv growth (1) or fixed Labor Pdtv (2))
b0=2 
#Type of Model (Labor Constrained (1) or Capital Constrained (2))
b1=2 

#PARAMETERS
b2=1 #a0
b3=0.025 #a_Gr
if(b0==1)
{b3=0.01} #Calibration b0=1}
b4=0.1 #a_Resp_K_Gr
b5=1 #a_Resp_Emp
b6=0.95 #a_Resp_Emp_Normal	
b8=0.9 #Consumption propensity
b10=0.5 #Mark_Up_Rt(0)
b12=1 #Inv_Resp_To_K_Prod_Grwth
if(b0==1)
{b12=0} #Calibration b0=1
b13=0.06 #Inv_Resp_To_K_Util
b14=3 #K/Q
b15=0.8 #Target Capacity Util Rt
b16=0.1 #dep_rt
if(b0==1)
{b16=0.05} #Calibration b0=1
b18=400 #K0
b19=100 #L
b21=0.1 #MarkUp_Adjstmt_Parameter
b22=1 #Target_Demand/Supply_Rt	1
b24=0 #HH_Autonomous_Resp_to_Wealth	0
b26=0 #Gov_Autonomous_Spending_Rt	0

# Declareing vectors
t0=1
I=Ir=u=K=D=Y=Ku=sigma=Sc=Sl=M=Khat=a=ahat=e=vector()

#Initial conditions
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
a[t0]=b2 #Labor Productivity
Khat[t0]=0 #Capital stock growth
ahat[t0]=b2 #Labor productivity growth
e[t0]=Y[t0]/(b19*a[t0]) #Total Employment

# Dynamic Equations  
for(t in 2:1000) {
  K[t]=K[t-1]*(1-b16+Ir[t-1])  
  u[t]=Ku[t-1]/K[t-1]
  
  if (t!=2) {
    if (b1==2) 
      {M[t]=M[t-1]*(1+b21*(D[t-1]/Sc[t-1]-b22))} 
    else 
      {M[t]=M[t-1]*(1+b21*(D[t-1]/Sl[t-1]-b22))}
  }
  
  Khat[t]=(K[t]-K[t-1])/K[t-1]
  a[t]=a[t-1]*(1+b3+b4*Khat[t-1]+b5*(Y[t-1]/(a[t-1]*b19)-b6))
  ahat[t]=b3+b4*Khat[t-1]+b5*(Y[t-1]/(a[t-1]*b19)-b6)
  
  if (b0==2)  
    {Ir[t]=b16+b13*(u[t]-b15)+b12*b3}
  else 
  {Ir[t]=b16+b13*(u[t]-b15)+b12*ahat[t-1]}

  
  I[t]=Ir[t]*K[t]
  Sc[t]=K[t]/b14
  Sl[t]=a[t]*b19  
  sigma[t]=1/(1+M[t])
  D[t]=I[t]/(1-b8*sigma[t])
  Y[t]=D[t]
  Ku[t]=b14*Y[t]
  e[t]=Y[t]/(b19*a[t])
}

#Graphs
par(mfrow=c(1,2), mar=c(3.1, 2.9, 2.5, 1), mgp=c(2, 1, 0), las=0, cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
plot( K, type="l", main="Total Capital", xlab="t", ylab="K",col="blue" )
plot( Khat, type="l", main="Capital Growth Rate", xlab="t", ylab="Khat",col="blue" )
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
plot( ahat, type="l", main="Productivity Growth", xlab="t", ylab="ahat",col="blue" )
plot( e, type="l", main="Employment rate", xlab="t", ylab="e",col="blue" )

# Printing the results in a .pdf file
if (print_results==1)
{
pdf("Results_Onder_Model.pdf") 
  # Graph Characteristics
  # par(mfrow=c(1,1), mar=c(3.1, 2.9, 2.5, 1), mgp=c(2, 1, 0), las=0, cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), lwd = 2, mgp=c(3, 1, 0), las=0, cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
  
  #Graphs
  plot( K, type="l", main="Total Capital", xlab="t", ylab="K",col="blue" )
  plot( Khat, type="l", main="Capital Growth Rate", xlab="t", ylab="Khat",col="blue" )
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
  #plot( a, type="l", main="Productivity", xlab="t", ylab="a",col="blue" )
  plot( ahat, type="l", main="Productivity Growth", xlab="t", ylab="ahat",col="blue" )
  plot( e, type="l", main="Employment rate", xlab="t", ylab="e",col="blue" )
#Close the pdf file
dev.off() 
}




############################### DIFFERENTIAL EQUATIONS ###########################3


OnderModel <- function (t, y, parms) {
  with(as.list(y), {
    
    #Differential Equations
    dK=K*(-b16+g)  
    dM= M*(b21*((((g*K)/(1-b8*(1/(1+M))))/(K/b14))-b22))
    dg= b16+b13*((((g*K)/(1-b8*(1/(1+M))))/(K/b14))-b15)+b12*b3-g
    da=b3+b4*K+b5*((g*K)/(1-b8*(1/(1+M)))/(a*b19)-b6)
    list(c(dK, dM, dg, da)) })

}

#Initial Conditions
yini <- c(K=b18, M=b10, g=b16, a=b2)

#Computing the dynamic equation (ODE)
times <- seq(from = 0, to = 100, by = 0.1)
out <- ode(y = yini, times = times, func =OnderModel ,
           parms = NULL)
par(mar=c(2, 2, 2, 2))
par(mfrow=c(2, 2))

#plot(out, lwd = 2, main="title")

#Correctly graph the equations
Kde=out[,2]
Mde=out[,3]
gde=out[,4]
ade=out[,5]

##

Kdehat=(Kde[-1]-Kde[-length(Kde)])/Kde[-length(Kde)]
Yde=Dde=Kde*gde/(1-b8*(1/(1+Mde)))
Ydehat=(Yde[-1]-Yde[-length(Yde)])/Yde[-length(Yde)]
adehat=(ade[-1]-ade[-length(Yde)])/ade[-length(Yde)]
edehat=Ydehat-adehat

##

plot(Kde, lwd = 2, type='l', main="Capital Stock", xlab="time", ylab="K")
plot(Mde, lwd = 2, type='l', main="Mark-Up", xlab="time", ylab="M")
plot(gde, lwd = 2, type='l', main="Investment rate", xlab="time", ylab="Ir")
plot(adehat, lwd = 2, type='l', main="Productivity growth", xlab="time", ylab="ahat")
plot(edehat, lwd = 2, type='l', main="Employment growth", xlab="time", ylab="ehat")
plot(Kdehat, lwd = 2, type='l', main="Capital Stock Growth", xlab="time", ylab="Khat")
plot(Ydehat, lwd = 2, type='l', main="Output Growth", xlab="time", ylab="Yhat")

