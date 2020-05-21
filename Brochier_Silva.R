############## Labor Market ################# 
########### Brochier & Silva (2019) ######### 
###### Cambridge Journal of Economics #######
#############################################

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
tau=0.37 #Ratio of taxes on household income
sigma=0.34 #Share of government expenditures on income
alfa1=0.8 #Fraction of after-tax wages - households consumption
alfa2=0.03375 #Fraction of stock of wealth - households consumption
a=0.1 #Equities are a fixed proportion of the capital stock at the beginning of the period
x=0.001 #capacity utilization band
lambda0=0.08 #proportion of household wealth allocated in equities
pe=0.9880292 #Prices of Equities
delta=0.044 #Depreciation rate
ir=0.02 #Interest rate
upsilon=2.5 #Capital-output ratio
mu=0.7 #Mark-up
pi=0.411765 # Profit share
#pi=mu/(1+mu)
sf=0.4 #Firms retain a fraction of their profit
gamma=0.014 #Speed of adjustment of the propensity to invest to the discrepancies between the actual utilization rate and the desired utilization rate.
h=0.2 #Investment rate
un=0.8 #Natural capacity utilization rate
lambda=lambda0-ir #proportion of household wealth allocated in equities

#Long-run Variables
#m_star=l_star=0.7953325
#b_star=0.7525772
#upsilonh_star=1.646715
#g_star=
#K0=100

# Declaring vectors
t0=1
B=G=Y=Te=Yh=W=FD=L=C=I=Vh=Yd=Sh=FU=M=K=E=Yfc=upsilonh=u=gk=h=gy=u2=pe=Fe=Fg=rn=rg=gvh=vector()
nd=ns=n=beta=betahat=e=ge=w=omegaw=omegaf=omega=phat=gu=What=p=phat=gu=vector()


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
M[t0]=10 #Deposit share in wealth
K[t0]=250 #Capital Stock
E[t0]=a*K[t0] #Equities
Yfc[t0]=K[t0]/upsilon  #Full capacity output 
upsilonh[t0]=upsilon #normalized household wealth to capital ratio
u[t0]=Y[t0]/Yfc[t0] #Capacity utilization rate
h[t0]=0.2 #Investment rate
pe[t0]=lambda*Vh[t0]/E[t0] #Price of equities
gk[t0]=0.2 #Growth of capital stock
gy[t0]=gk[t0] #Growth of output
Fg[t0]=pi*Y[t0] #Gross profit
rg[t0]=pi*u[t0]/upsilon #Net profit rate
gvh[t0]=0.02 #Growth of wealth

# Dynamic Equations  
for(t in 2:1000) 
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

  Fe[t]=pi*Y[t]-L[t-1]*ir
  Fg[t]=pi*Y[t]
  #rn[t]=pi*(u[t]/upsilon)-ir*l[t-1]/(1+gk[t-1])
  rg[t]=pi*u[t]/upsilon

  ############## Labor Market ################# 
  #####Brochier(2020) in Metroeconomica #######
  #############################################
  eta0=0.01
  eta1=0.01
  xi=0.01
  epsilon1=0.1
  epsilon2=0.1
  epsilon3=0.1
  phi1=0.1
  phi2=0.1
  psi=0.1
  theta0=0
  theta1=0.001
  
  beta[t0]=1
  betahat[t0]=0.02
  nd[t0]=Y[t0]/beta[t0]
  nd[t0]=80
  ns[t0]=100
  omegaw[t0]=0.7
  omegaf[t0]=0.6
  ge[t0]=0.02
  w[t0]=60 #Wages
  What[t0]=0.01
  phat[t0]=0.01
  p[t0]=100
  omega[t0]=w[t0]/(beta[t0]*p[t0])
  e[t0]=0.8  
  gn[t0]=theta0+theta1*e[t0]

  beta[t]=beta[t-1]*(1+betahat[t-1])
  betahat[t]=eta0+eta1*gy[t-1]
 
  w[t]=w[t-1]*(1+What[t-1])
  p[t]=p[t-1]*(1+phat[t-1])
  
  omega[t]=w[t]/(beta[t]*p[t])
  
  e[t]=nd[t-1]/ns[t-1]  # emplyment explodes
  ge[t]=(e[t]-e[t-1])/e[t-1]
  omegaw[t]=omegaw[t-1]+xi*ge[t]
  What[t]=epsilon1*(omegaw[t]-omega[t])+epsilon2*phat[t-1]+epsilon3*betahat[t]
  
  gu[t]=(u[t]-u[t-1])/u[t-1]
  omegaf[t]=omegaf[t-1]-psi*gu[t]
  phat[t]=phi1*(omega[t]-omegaf[t])+phi2*(What[t-1])
  
  gn[t]=theta0+theta1*e[t-1]
  nd[t]=Y[t]/beta[t]
  ns[t]=ns[t-1]*(1+gn[t])
    
  ######Dynamic simultaneous system#####
  
  #gy[t]=(Y[t]-Y[t-1])/Y[t-1]
  #gk[t]=(h[t]*u[t]/upsilon)-delta
  #h[t]=h[t-1]+h[t-1]*gamma*(u[t]-un)
  #u[t]=(1/((1+gk[t-1])*(1-h[t]-alfa1*(1-tau)*(1-pi)-(sigma/(1+gy[t-1])))))*alfa2*upsilonh[t-1]*upsilon
  
  #####################################
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

