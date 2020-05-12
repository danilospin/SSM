library(deSolve)
library(rootSolve)
library(psych)
library(scatterplot3d)

#List of Parameters of the Model:
mu=0.8 #desired utilization rate
gamma=0.1 #adjustment speed of investment wrt u
gz=0.02 #growth of exogenous spending
s=0.45 # propensity to save
v=2 #(desired)capital-output ratio
delta=0.05 #depreciation rate of capital
rho = 0.02 #labour productivity growth rate
xi = 0.2 # endogenous part of productivity growth
ebar = 0.8
lambda = 0.2
zeta=0.0

#h is the fraction of income that is invested
#u is the capital utilization rate

ssmmod <- function (t, y, parms) {
  with(as.list(y), {
    
    #Differential Equations
    dh=h*gamma*(u-mu)
    du=u*(gz+(h*gamma*(u-mu))/(s-h)-u*h/v+delta)
    de=e*(gz+(h*gamma*(u-mu))/(s-h)-r)
    
    
    #Bart's Version
    #dr=(xi*(u*h*gamma*(u-mu)+h*u*(gz+(h*gamma*(u-mu))/(s-h)-u*h/v+delta))/v)+lambda*(e-ebar)
    #Onder's Version'
    dr=r*(zeta+xi*(u*h)/v-delta)+lambda*(e-ebar)
    list(c(dh, du, de, dr)) })
}

#Initial Conditions
yini <- c(h = 0.2, u=0.8, e=0.7, r=0.02)

#Computing the dynamic equation (ODE)
times <- seq(from = 0, to = 500, by = 0.1)
out <- ode(y = yini, times = times, func = ssmmod,
           parms = NULL)
par(mar=c(2, 2, 2, 2))
par(mfrow=c(2, 2))
#plot(out, lwd = 2, main="title")

#Correctly graph the equations
hvar=out[,2]
uvar=out[,3]
evar=out[,4]
rvar=out[,5]

gk = hvar*uvar/v-delta
g = gz+(hvar*gamma*(uvar-mu))/(s-hvar)
mult = 1/(s-hvar)

plot(hvar, lwd = 2, type='l', main="Investment ratio", xlab="time", ylab="h")
plot(uvar, lwd = 2, type='l', main="Utilization rate", xlab="time", ylab="u")
plot(evar, lwd = 2, type='l', main="Employment rate", xlab="time", ylab="e")
plot(rvar, lwd = 2, type='l', main="Productivity growth", xlab="time", ylab="rho")

plot(gk, lwd = 2, type='l', main="Growth of capacity", xlab="time", ylab="gK")
plot(g, lwd = 2, type='l', main="Growth of output", xlab="time", ylab="g")
plot(mult, lwd = 2, type='l', main="Multiplier", xlab="time", ylab="m")

stab = 1-s+v*(gz+delta)/mu+gamma*v
h_eq = v*(gz+delta)/mu

res=cbind(out[,1], gk, g, mult)
setwd("C:/Users/Danilo/Desktop/Paper Supermultiplier")
write.csv(out, "results.csv")
write.csv(res, "results2.csv")

#Jacobian
#jacobiansteady=runsteady(y=yini, time = c(0, Inf), func=ssmmod, parms=NULL)
#Steadystate=jacobiansteady$y
#jacobiano=jacobian.full(y = jacobiansteady$y, func = ssmmod, parms=NULL)
#Steadystate
#eigen(jacobiano)$values
