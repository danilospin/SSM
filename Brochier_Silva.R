## Brochier-Silva ##

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


#PARAMETERS

# Declareing vectors
t0=1
a=vector()

#Parameters Baseline

tau=0.37
sigma=0.34
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
sf=0.4
gamma=0.014
h=0.2
un=0.8

#Long-run Variables
m_star=l_star=0.7953325
b_star=0.7525772
upsilonh_star=1.646715
g_star=
K0=100

  
#Initial conditions
B[t0]=
G[t0]=
Te[t0]=
Y[t0]=
W[t0]=
M[t0]=
vh[t0]=
FD[t0]=
Yh[t0]=
  Yd[t0]=
  C[t0]=
sh[t0]=
pe[t0]= 
  u[t0]=
  Yfc[t0]=
  gk[t0]=
  h[t0]=
  
# Dynamic Equations  
for(t in 2:1000) {
  B[t]=B[t-1]+G[t]-Te[t]+ir*B[t-1]
  G[t]=sigma*Y[t-1]
  Te[t]=tau*Yh[t]
  Yh[t]=W[t]+FD[t]+ir*(B[t-1]+M[t-1])
  W[t]=(1-pi)*Y[t]
  Yd[t]=(1-tau)*Yh[t]
  C[t]=alfa1*(1-tau)*W[t]+alfa2*Vh[t-1]
  Sh[t]=Yd[t]-C[t]
  
  lambda=lambda0-ir
  pe[t]=lambda*Vh[t]/E[t]
  
  Vh[t]=Vh[t-1]+Sh+(pe[t]-pe[t-1])*E[t-1]
  M[t]=M[t-1]+Sh-pe*(E[t]-E[t-1])
  
  pi=mu/(1+mu)
  I[t]=h*Y[t]
  h[t]=h[t-1]+h[t-1]*gamma*(u[t]-un)
  
  K[t]=K[t-1]-delta*K[t-1]+I[t]
  Yfc=K[t-1]/upsilon
  
  u[t]=Y[t]/Yfc[t]
  gk[t]=(h[t]*u[t]/upsilon)-delta
  
  L[t]-L[t-1]+I[t]-FU[t]-pe[t]*(E[t]-E[t-1])
  E[t]=a*K[t-1]
  FU[t]=sf*(pi*Y[t]-ir*L[t-1])
  FD[t]=(1-sf)*(pi*Y[t]-ir*L[t-1])
  Fe[t]=pi*Y[t]-ir*L[t-1]
  Fg[t]=pi*Y[t]
  
  rn[t]=pi*(u[t]/upsilon)-ir*l[t-1]/(1+gk[t-1])
  rg[t]=pi*u/upsilon
  
  Y[t]=C[t]+I[t]+G[t]
  
  upsilonh[t]=Vh[t]/K[t-1]
  u[t]=(1/(1+gk[t-1])*(1-h-alfa1*(1-tau)*(1-pi)-(sigma/(1+gy(t-1)))))*alfa2*upsilonh(t-1)*upsilon
}


#Graphs
#par(mfrow=c(1,2), mar=c(3.1, 2.9, 2.5, 1), mgp=c(2, 1, 0), las=0, cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
#plot( K, type="l", main="Total Capital", xlab="t", ylab="K",col="blue" )
#plot( Khat, type="l", main="Capital Growth Rate", xlab="t", ylab="Khat",col="blue" )
#plot( Ku, type="l", main="Used Capital", xlab="t", ylab="Ku",col="blue" )

# Printing the results in a .pdf file
if (print_results==1)
{
  pdf("Results_Brochier_Silva.pdf") 
  # Graph Characteristics
  # par(mfrow=c(1,1), mar=c(3.1, 2.9, 2.5, 1), mgp=c(2, 1, 0), las=0, cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), lwd = 2, mgp=c(3, 1, 0), las=0, cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
  
  #Graphs
  #plot( K, type="l", main="Total Capital", xlab="t", ylab="K",col="blue" )
  plot( Khat, type="l", main="Capital Growth Rate", xlab="t", ylab="Khat",col="blue" )
  #plot( Ku, type="l", main="Used Capital", xlab="t", ylab="Ku",col="blue" )
  plot( u, type="l", main="Utilization Capacity", xlab="t", ylab="u",col="blue" )
  plot( M, type="l", main="Mark-up", xlab="t", ylab="M",col="blue" )
  plot( sigma, type="l", main="w/a", xlab="t", ylab="w/a",col="blue" )
  plot( Ir, type="l", main="Investment Rate", xlab="t", ylab="Ir",col="blue" )
  #plot( I, type="l", main="Investment", xlab="t", ylab="I",col="blue" )
  #plot( Sc, type="l", main="A. Supply (Capital)", xlab="t", ylab="Sc",col="blue" )
  #plot( Sl, type="l", main="A. Supply (Labor)", xlab="t", ylab="Sl",col="blue" )
  #plot( D, type="l", main="Aggregate Demand", xlab="t", ylab="D",col="blue" )
  plot( D/Sc, type="l", main="Demand/C Supply Ratio", xlab="t", ylab="D/S",col="blue" )
  plot( D/Sl, type="l", main="Demand/L Supply Ratio", xlab="t", ylab="D/S",col="blue" )
  #plot( a, type="l", main="Productivity", xlab="t", ylab="a",col="blue" )
  plot( ahat, type="l", main="Productivity Growth", xlab="t", ylab="ahat",col="blue" )
  plot( e, type="l", main="Employment rate", xlab="t", ylab="e",col="blue" )
  #Close the pdf file
  dev.off() 
}




############################### DIFFERENTIAL EQUATIONS ###########################3


BrochierModel <- function (t, y, parms) {
  with(as.list(y), {
    
    #Differential Equations
    dK=K*(b16+((b16+b12*b3)*(1-b13*(b14*(K)/(1-b8*(1/(1+M)))/K-b15))^(-1)))
    dM=M*(b21*(((((b16+b12*b3)* (1-b13*(b14*(K)/(1-b8*(1/(1+M)))/K-b15))^(-1))*K)/(1-b8*(1/(1+M))))/(K/b14)-b22))
    
    
    list(c(dK, dM)) })
}

#Initial Conditions
yini <- c(K= 300, M=1)

#Computing the dynamic equation (ODE)
times <- seq(from = 0, to = 1000, by = 0.1)
out <- ode(y = yini, times = times, func =BrochierModel ,
           parms = NULL)
par(mar=c(2, 2, 2, 2))
par(mfrow=c(2, 2))
plot(out, lwd = 2, main="title")

