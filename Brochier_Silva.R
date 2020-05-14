## Brochier-Silva ##

#Setting the libraries
library(deSolve)
library(rootSolve)
library(psych)
library(scatterplot3d)

#Setting the Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

#SETTING THE MODEL
#print results in a .pdf file? 1 - Yes, 2 - No
print_results=1

#PARAMETERS
tau=0.37
sigma=0.34 #Share of government expenditures on income
alfa1=0.8
alfa2=0.03375
a=0.1
x=0.001
lambda0=0.08
pe=0.9880292
delta=0.044
ir=0.02
upsilon=2.5
mu=0.7
pi=0.411765
#pi=mu/(1+mu)
sf=0.4
gamma=0.014
h=0.2
un=0.8
lambda=lambda0-ir

#Long-run Variables
#m_star=l_star=0.7953325
#b_star=0.7525772
#upsilonh_star=1.646715
#g_star=
#K0=100

# Declaring vectors
t0=1
B=G=Y=Te=Yh=W=FD=L=C=I=Vh=Yd=Sh=FU=M=K=E=Yfc=upsilonh=u=gk=h=gy=u2=pe=Fe=Fg=rn=rg=gvh=vector()

#Initial conditions
L[t0]=50 #Financing Investment
FD[t0]=0 #Financial income to households (interest on deposits and bills and dividends)
Yh[t0]=60 #Household income
Te[t0]=tau*Yh[t0] #Taxation of household income
Y[t0]=80 #Aggregate income
W[t0]=(1-pi)*Y[t0] #Wages
C[t0]=alfa1*(1-tau)*W[t0] #Total Consumption
B[t0]=0 #Government debt
G[t0]=sigma*Y[t0] #Government expenditures
Vh[t0]=100 #Stock of wealth
Yd[t0]=(1-tau)*Yh[t0] #Disposable income
Sh[t0]=Yd[t0]-C[t0] #Household savings
FU[t0]=0 #Firms retain a fraction of their profit discounting the payment of interest on loans
M[t0]=10 #Deposits
K[t0]=250 #Capital Stock
E[t0]=a*K[t0] #Equities
Yfc[t0]=K[t0]/upsilon  #Full capacity output 
upsilonh[t0]=upsilon #normalized household wealth to capital ratio
u[t0]=Y[t0]/Yfc[t0]
u2[t0]=u[t0]
h[t0]=0.2
pe[t0]=lambda*Vh[t0]/E[t0]
gk[t0]=0.2
gy[t0]=gk[t0]
Fg[t0]=pi*Y[t0]
rg[t0]=pi*u[t0]/upsilon
gvh[t0]=0.02

# Dynamic Equations  
for(t in 2:100) 
  {
  #u[t]=(1/((1+gk[t-1])*(1-h[t]-alfa1*(1-tau)*(1-pi)-(sigma/(1+gy[t-1])))))*alfa2*upsilonh[t-1]*upsilon
  #h[t]=h[t-1]+h[t-1]*gamma*(u[t]-un)
  
#I used the LAG to avoid overdetermination - but this is wrong
  W[t]=(1-pi)*Y[t-1] #used Lag 
  I[t]=h*Y[t-1] #Double determination, Investment uses lag of income?
  C[t]=alfa1*(1-tau)*W[t]+alfa2*Vh[t-1]
  G[t]=sigma*Y[t-1]
  Y[t]=G[t]+C[t]+I[t]
 #Y[t]=(1/(1-alfa1*(1-tau)*(1-pi)-h[t]))*(sigma*Y[t-1]+alfa2*Vh[t-1]) //Correct according to the paper
 #Y[t]=(1/(1-alfa1*(1-tau)*(1-pi)-h[t]-sigma))*(alfa2*Vh[t-1]) 
  
  E[t]=a*K[t-1]
  
  FU[t]=sf*(pi*Y[t]-ir*L[t-1])
  L[t]=L[t-1]+I[t]-FU[t]-(E[t]-E[t-1])#*pe[t]
  FD[t]=(1-sf)*(pi*Y[t]-ir*L[t-1])
  
  Yh[t]=W[t]+FD[t]+ir*(B[t-1]+M[t-1])
  Yd[t]=(1-tau)*Yh[t]
  Sh[t]=Yd[t]-C[t]
  
  K[t]=K[t-1]-delta*K[t-1]+I[t]
  #gk[t]=(K[t]-K[t-1])/K[t-1]
  
  #Vh[t]=Vh[t-1]+Sh[t]+E[t-1]#*(pe[t]-pe[t-1])
  
  Vh[t]=((1-(E[t-1]*lambda/E[t]))^(-1))*(Vh[t-1]+Sh[t]-pe[t-1]*E[t-1])
  gvh[t]=(Vh[t]-Vh[t-1])/Vh[t-1]
  pe[t]=lambda*Vh[t]/E[t]
  
  M[t]=M[t-1]+Sh[t]-(E[t]-E[t-1])*pe
  
  Te[t]=tau*Yh[t]
  B[t]=B[t-1]+G[t]-Te[t]+ir*B[t-1]
  
  upsilonh[t]=Vh[t]/K[t-1]
  
  Yfc[t]=K[t-1]/upsilon
  u[t]=Y[t]/Yfc[t]
  h[t]=h[t-1]+h[t-1]*gamma*(u[t]-un)
  gk[t]=(h[t]*u[t]/upsilon)-delta
  gy[t]=(Y[t]-Y[t-1])/Y[t-1]

  u2[t]=(1/((1+gk[t-1])*(1-h[t]-alfa1*(1-tau)*(1-pi)-(sigma/(1+gy[t-1])))))*alfa2*upsilonh[t-1]*upsilon

  Fe[t]=pi*Y[t]-L[t-1]*ir
  Fg[t]=pi*Y[t]
  #rn[t]=pi*(u[t]/upsilon)-ir*l[t-1]/(1+gk[t-1])
  rg[t]=pi*u[t]/upsilon
  
}

#Graphs
par(mfrow=c(1,2), mar=c(3.1, 2.9, 2.5, 1), mgp=c(2, 1, 0), las=0, cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
plot( gy, type="l", main="Growth", xlab="t", ylab="gy",col="blue" )
plot( u, type="l", main="Capacity", xlab="t", ylab="u",col="blue" )
plot( gvh, type="l", main="Growth Wealth", xlab="t", ylab="gvh",col="blue" )
#plot( Khat, type="l", main="Capital Growth Rate", xlab="t", ylab="Khat",col="blue" )
#plot( Ku, type="l", main="Used Capital", xlab="t", ylab="Ku",col="blue" )

# Printing the results in a .pdf file
#if (print_results==1)
#{
#  pdf("Results_Brochier_Silva.pdf") 
# Graph Characteristics
# par(mfrow=c(1,1), mar=c(3.1, 2.9, 2.5, 1), mgp=c(2, 1, 0), las=0, cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
#  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), lwd = 2, mgp=c(3, 1, 0), las=0, cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
  #Graphs
  #plot( K, type="l", main="Total Capital", xlab="t", ylab="K",col="blue" )
  #plot( Khat, type="l", main="Capital Growth Rate", xlab="t", ylab="Khat",col="blue" )
  #plot( Ku, type="l", main="Used Capital", xlab="t", ylab="Ku",col="blue" )
  #plot( u, type="l", main="Utilization Capacity", xlab="t", ylab="u",col="blue" )
  #plot( M, type="l", main="Mark-up", xlab="t", ylab="M",col="blue" )
  #plot( sigma, type="l", main="w/a", xlab="t", ylab="w/a",col="blue" )
  #plot( Ir, type="l", main="Investment Rate", xlab="t", ylab="Ir",col="blue" )
  #plot( I, type="l", main="Investment", xlab="t", ylab="I",col="blue" )
  #plot( Sc, type="l", main="A. Supply (Capital)", xlab="t", ylab="Sc",col="blue" )
  #plot( Sl, type="l", main="A. Supply (Labor)", xlab="t", ylab="Sl",col="blue" )
  #plot( D, type="l", main="Aggregate Demand", xlab="t", ylab="D",col="blue" )
  #plot( D/Sc, type="l", main="Demand/C Supply Ratio", xlab="t", ylab="D/S",col="blue" )
  #plot( D/Sl, type="l", main="Demand/L Supply Ratio", xlab="t", ylab="D/S",col="blue" )
  #plot( a, type="l", main="Productivity", xlab="t", ylab="a",col="blue" )
  #plot( ahat, type="l", main="Productivity Growth", xlab="t", ylab="ahat",col="blue" )
  #plot( e, type="l", main="Employment rate", xlab="t", ylab="e",col="blue" )
  #Close the pdf file
  #dev.off() 
#}

############################### DIFFERENTIAL EQUATIONS ###########################3

#BrochierModel <- function (t, y, parms) {
#  with(as.list(y), {
#    
#Differential Equations
#    dK=K*(b16+((b16+b12*b3)*(1-b13*(b14*(K)/(1-b8*(1/(1+M)))/K-b15))^(-1)))
#    dM=M*(b21*(((((b16+b12*b3)* (1-b13*(b14*(K)/(1-b8*(1/(1+M)))/K-b15))^(-1))*K)/(1-b8*(1/(1+M))))/(K/b14)-b22))
    
#    list(c(dK, dM)) })
#}

#Initial Conditions
#yini <- c(K= 300, M=1)

#Computing the dynamic equation (ODE)
#times <- seq(from = 0, to = 1000, by = 0.1)
#out <- ode(y = yini, times = times, func =BrochierModel ,
#           parms = NULL)
#par(mar=c(2, 2, 2, 2))
#par(mfrow=c(2, 2))
#plot(out, lwd = 2, main="title")

