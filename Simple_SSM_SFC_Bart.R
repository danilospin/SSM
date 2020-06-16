library(deSolve)
library(rootSolve)
library(psych)
library(scatterplot3d)

#List of Parameters of the Model:
mu=0.9 #desired utilization rate
gamma=0.01 #adjustment speed of investment wrt u
sw = 0.175 # propensity to save, workers
v=2 #(desired)capital-output ratio
gre = 0.02 # exogenous productivity growth
delta = 0.05 # depreciation
zeta = 0.05 # share of W consumed
lambda = 1.0 #responsiveness of sigma (s)
ebar = 0.8 #neutral employment rate

#h is the fraction of income that is invested
#u is the capital utilization rate

ssmmod <- function (t, y, parms) {
  with(as.list(y), {
    
    #Differential Equations
    dh=h*gamma*(u-mu)
    du=u*((-zeta*(1-s-h)/(1-s+s*sw-h))+h*gamma*(u-mu)/(1-s+s*sw-h)-u*h/v+delta)
    dw=-w*zeta*(1-s-h)/(1-s+s*sw-h)
    ds=-s*(lambda*(e-ebar)+gre)
    de=e*(-zeta*(1-s-h)/(1-s+s*sw-h)+(-s*(lambda*(e-ebar)+gre)*(1-sw)+h*gamma*(u-mu))/(1-s+s*sw-h)-gre)
    list(c(dh, du, dw, ds, de)) })
}


eq_h = (delta+gre)*v/mu
eq_u = mu
eq_s = (zeta+gre)*(1-(delta+gre)*v/mu)/(zeta+gre*(1-sw))
eq_e = ebar - gre/lambda
eq_g = -zeta*(1-eq_s-eq_h)/(1-eq_s+eq_s*sw-eq_h)


#Initial Conditions
yini <- c(h=0.05, u=0.85, w=1, s=0.75, e=0.84)

#Computing the dynamic equation (ODE)
times <- seq(from = 0, to = 500, by = 1)
out <- ode(y = yini, times = times, func = ssmmod,
           parms = NULL)
plot(out, lwd = 2, main="title")

#Correctly graph the equations
hvar=out[,2]
uvar=out[,3]
wvar=out[,4]
svar=out[,5]
evar=out[,6]

mult = 1/(1-svar+svar*sw-hvar)

par(mfrow=c(1,2), mar=c(3.1, 2.9, 2.5, 1), mgp=c(2, 1, 0), las=0, cex.lab=0.7, cex.axis=0.7, cex.main=0.7, cex.sub=0.7)

plot(uvar, lwd = 2, type='l', main="Utilization rate", xlab="time", ylab="u")
plot(svar, lwd = 2, type='l', main="Workers income share", xlab="time", ylab="sigma")
plot(hvar, lwd = 2, type='l', main="Investment share", xlab="time", ylab="h")
plot(evar, lwd = 2, type='l', main="Employment rate", xlab="time", ylab="e")

plot(wvar, lwd = 2, type='l', main="Workers assets", xlab="time", ylab="w")
plot(mult, lwd = 2, type='l', main="Multiplier", xlab="time", ylab="m")

#Jacobian
#jacobiansteady=runsteady(y=yini, time = c(0, Inf), func=ssmmod, parms=NULL)
#Steadystate=jacobiansteady$y
#jacobiano=jacobian.full(y = jacobiansteady$y, func = ssmmod, parms=NULL)
#Steadystate
#eigen(jacobiano)$values
