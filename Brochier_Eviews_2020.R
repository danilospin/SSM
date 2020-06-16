############## Labor Market ################# 
############ Brochier (2020) ################
############# METROECONOMICA ################
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
tau=0.35 #Ratio of taxes on household income
sigma=0.34 #Share of government expenditures on income
alpha1=0.8 #Fraction of after-tax wages - households consumption
alpha2=0.065613 #Fraction of stock of wealth - households consumption
zeta=0.1 #Equities are a fixed proportion of the capital stock at the beginning of the period
x=0.0001 #capacity utilization band
lambda=0.44
lambda0=0.449901 #proportion of household wealth allocated in equities
delta=0.077 #Depreciation rate
i=0.0200109
ir=0.00990099 #Interest rate
v=1.5
pi=0.4444445 # Profit share
sf=0.3 #Firms retain a fraction of their profit
gamma=0.05 #Speed of adjustment of the propensity to invest to the discrepancies between the actual utilization rate and the desired utilization rate.
h=0.18 #Investment rate
un=0.8499798 #Natural capacity utilization rate
theta0=-0.065
theta1=0.1
eta0=0
eta1=0
xi=0
psi=0
epsilon1=0.5
epsilon2=0.1
epsilon3=0
phi1=0.5
phi2=0.1
sf=0.3
lambda=lambda0-ir #C18

# Declaring vectors
t0=1
B=G=Y=Te=Yh=W=FD=L=C=I=Vh=Yd=Sh=FU=M=K=Eq=Yfc=upsilonh=u=gk=h=gy=u2=pe=Fe=Fg=sh=rn=rg=gvh=y=WB=yh=vh=gns=vector()
nd=ns=n1=n=beta=betahat=e=ge=w=omegaw=omegaf=omega=phat=gu=Fn=What=p=phat=gu=k=yd=te=yfc=id=gd=c=UC=vector()

#Initial conditions
p[t0]=1
y[t0]=80
Y[t0]=80 #Aggregate income
W[t0]=(1-pi)*Y[t0] #Wages
beta[t0]=0.01
nd[t0]=0.8

WB[t0]=n[t0]*W[t0]
C[t0]=alpha1*(1-tau)*W[t0] #Total Consumption
c[t0]=C[t0]/p[t0] # added
yfc[t0]=100 # added
u[t0]=y[t0]/yfc[t0] #Capacity utilization rate
h[t0]=0.2 #Investment rate
id[t0]=h[t0]*y[t0] # added
I[t0]=30
gd[t0]=10
G[t0]=gd[t0]*p[t0] #Government expenditures
k[t0]=0.01
e[t0]=0.8
ge[t0]=0.01 # added
omegaw[t0]=0.573575
gu[t0]=0.01 # added
omegaf[t0]=0.537536
gy[t0]=0.01 # added
betahat[t0]=0.02
gns[t0]=0.2
ns[t0]=0.82
omega[t0] =WB[t0]/Y[t0] # added
What[t0]=0.01
phat[t0]=0.01
te[t0] = 0.1 # added
Te[t0]=te[t0]*p[t0] #Taxation of household income
i[t0]= 0.01 # added
B[t0]=0 #Government debt
Fn[t0]=10 # added
FD[t0]=0 #Financial income to households (interest on deposits and bills and dividends)
Yh[t0]=60 #Household income
yh[t0]=Yh[t0]/p[t0]
Yd[t0]=(1-tau)*Yh[t0] #Disposable income
yd[t0] = Yd[t0]/p[t0] # added
Sh[t0]=Yd[t0]-C[t0] #Household savings
sh[t0] = yd[t0]-c[t0] # added
Eq[t0]=1 #Equities
Vh[t0]=100 #Stock of wealth
vh[t0]=Vh[t0]/p[t0]               
pe[t0]=lambda*Vh[t0]/Eq[t0] #Price of equities
K[t0]=250 #Capital Stock
FU[t0]=0.2 #Firms retain a fraction of their profit discounting the payment of interest on loans
M[t0]=10 #Deposit share in wealth
L[t0]=50 #Financing Investment
Fe[t0]=Y[t0]-WB[t0]
UC[t0]=W[t0]/beta[t0]

gy[t0]=0.02 #Growth of output
i[t0]= 0.01 # added
n1[t0]=nd[t0]-ns[t0]
n[t0]=nd[t0]*(1-n1[t0])+ns[t0]*n1[t0]
#WB[t0]=n[t0]*W[t0] #added


for(t in 2:100) 
{
  p[t]=p[t-1]*(1+phat[t-1]) #C44
  y[t]=(alpha2*vh[t-1]+(y[t-1]*sigma))*(1-alpha1*(1-tau)*(1/beta[t-1]*(1+betahat[t-1]))*(W[t-1]*(1+What[t-1]))/p[t]-(h[t-1]+h[t-1]*gamma*(u[t-1]-un)))^(-1) ##PROBLEM
  Y[t]=y[t]*p[t]
  W[t]=W[t-1]*(1+What[t-1]) #C41
  beta[t]=beta[t-1]*(1+betahat[t-1]) #C36
  nd[t]=y[t]/beta[t] #C33
  gns[t]=theta0+theta1*e[t-1] #C35
  ns[t]=ns[t-1]*(1+gns[t]) #implied, but should it be gns[t]?
  n1[t]=nd[t]-ns[t]
  n[t]=nd[t]*(1-n1[t])+ns[t]*n1[t] #C33
  WB[t]=n[t]*W[t] #C38
  c[t]=alpha1*(1-tau)*WB[t]/p[t]+alpha2*vh[t-1] #C13
  C[t]=c[t]*p[t] #C14
  yfc[t]=k[t-1]/v #C25
  u[t]=y[t]/yfc[t] #C24
  h[t]=h[t-1]+h[t-1]*gamma*(u[t]-un)#*(1-dx) #C26
  id[t]=h[t]*y[t] #C21
  I[t]=id[t]*p[t] #C20
  gd[t]=y[t-1]*sigma #C1
  G[t]=gd[t]*p[t] #C3
  k[t]=k[t-1]+id[t]-delta*k[t-1] #C23
  e[t]=(nd[t-1]/ns[t-1])*(1-n1[t])+n1[t] #C34
  ge[t]=(e[t]-e[t-1])/e[t-1] #growth rate
  omegaw[t]=omegaw[t-1]+xi*ge[t] #C43
  gu[t]=(u[t]-u[t-1])/u[t-1] #growth rate
  omegaf[t]=omegaf[t-1]-psi*gu[t] #C46
  gy[t]=(y[t]-y[t-1])/y[t-1] #growth rate
  betahat[t]=eta0+eta1*gy[t-1] #C37
  omega[t]=WB[t]/Y[t] #C39
  What[t]=epsilon1*(omegaw[t]-omega[t])+epsilon2*phat[t-1]+epsilon3*betahat[t] #C42
  phat[t]=phi1*(omega[t]-omegaf[t])+phi2*(What[t-1]) #C45
  te[t]=tau*yh[t-1] #C2
  Te[t]=te[t]*p[t] #C4
  i[t]=(1+ir)*(1+phat[t])-1 #C6
  #ir=i[t]*p[t]
  #lambda[t]=lambda0*ir[t]
  B[t]=B[t-1]+G[t]-Te[t]+i[t]*B[t-1] #C5
  Fn[t]=Y[t]-WB[t]-i[t]*L[t-1] #C29
  FD[t]=(1-sf)*Fn[t] #C30
  Yh[t]=WB[t]+FD[t]+i[t]*(B[t-1]+M[t-1]) #C7
  yh[t]=Yh[t]/p[t] #C8
  Yd[t]=(1-tau)*Yh[t] #C9
  Sh[t]=Yd[t]-C[t] #C11
  sh[t]=Yd[t]/p[t]-c[t] #C12
  Eq[t]=zeta*k[t-1] #C32
  #pe[t]=lambda*Vh[t-1]/Eq[t] #It is Vh[t], but I Vh[t-1]
  ##Vh[t]=Vh[t-1]+Sh[t]+(pe[t]-pe[t-1])*Eq ##Doubt!
  Vh[t]=(Vh[t-1]+Sh[t]-pe[t-1]*Eq[t-1])*(1-lambda/Eq[t]*Eq[t-1])^(-1) #C15 PROBLEM!
  pe[t]=lambda*Vh[t]/Eq[t] #C19
  vh[t]=Vh[t]/p[t] #C16
  K[t]=k[t]*p[t] #C22
  FU[t]=sf*Fn[t] #C31
  M[t]=M[t-1]+Sh[t]-pe[t]*(Eq[t]-Eq[t-1])-(B[t]-B[t-1]) #C17
  L[t]=L[t-1]+I[t]-FU[t]-pe[t]*(Eq[t]-Eq[t-1]) #C27
  Fe[t]=Y[t]-WB[t] #C28
  yd[t]=Yd[t]/p[t] #C10
  UC[t]=W[t]/beta[t] #C40
}

#Graphs
par(mfrow=c(1,2), mar=c(3.1, 2.9, 2.5, 1), mgp=c(2, 1, 0), las=0, cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)
plot( gy, type="l", main="Growth", xlab="t", ylab="gy",col="blue" )
plot( gy, type="l", main="Growth", xlab="t", ylab="gy",col="blue" )
plot( u, type="l", main="Capacity", xlab="t", ylab="u",col="blue" )
plot( vh, type="l", main="Wealth", xlab="t", ylab="gvh",col="blue" )
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

