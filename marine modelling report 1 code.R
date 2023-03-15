#marine modelling rapport 

rm(list=ls())
library(deSolve)
library(spam)
library(viridis)
library(viridisLite)
library(fields)



param=c()
param$dz=1
param$z=seq(param$dz/2,200,by=param$dz)
n=length(param$z)
param$kp=0.1 #m^2/my mol N, self shading coeff.
param$kw=0.1 #m^-1 light absorbed by water.
param$I0=350  # mymol photons*m^2/sec
Hi=30 # maybe change the unit? mymol photons/m2/s halfsaturation of light
Hn=0.3# mymol/ m^3???
m=0.01*24#mortality, pr day
gmax=1.5 # pr dag
u=3# this is dropping of plankton parm. m/day
D=5*10^-5*3600*24 #m^2/day diffusivity.
nbot=30
alpha=1
eps=0.1 # day ^-1, remineralization rate.
bet=5 # strength of light effect, just a variable skalar.
gamma=0.2 # grazing m^3 (mmol N)^-1 day^-1

y=c(1:(n*3))*0
y[1]=1 # initial plankton, in nutrients.
y[(n+1):(2*n)]=20# initial nutrients
y[(2*n+1):(3*n)]=5# detritus.

lightemissions=function(P,Detrit,param,t){
  summerpleb=(cos(2*pi*(t/365))+1.01)*bet
# this is to set seasonal effect to not being there
  #summerpleb=1
  Q=param$kp*param$dz*((cumsum(P)-P/2)+(cumsum(Detrit)-Detrit/2))
  I=param$I0*exp(-param$kw*param$z-Q)*summerpleb
  return(I)
}




Ja=c() # advective flux
Jd=c() #diffusive flux
Jdn=c() # nutrient flux
Jdt=c()  #detritus flux
Jadt=c() # also detritus flux
 
P=c()
N=c()
dPdt=c()
dNdt=c()


diffusivefux=function(t,y,mu){
  #list(c(phi))
  P=y[1:n]
  N=y[(n+1):(n*2)]
  Detrit=y[(2*n+1):(3*n)]
  for (i in 2:n){
    Ja[i]=u*P[i-1]
    Jd[i]=-D*((P[i]-P[i-1])/param$dz)
    Jdn[i]=-D*((N[i]-N[i-1])/param$dz)
    Jadt[i]=u*Detrit[i-1]
    Jdt[i]=-D*(Detrit[i]-Detrit[i-1])/param$dz
    
  }
  Ja[c(1,n+1)]=0  
  Jd[c(1,n+1)]=0
  Jdn[c(1,n+1)]=c(0, -D*(nbot-N[n])/param$dz) 
  Jadt[c(1,n+1)]=c(0,u*Detrit[n])
  Jdt[c(1,n+1)]=0  
  J=Ja+Jd
  JD=Jadt+Jdt
  
  
  ii=lightemissions(P,Detrit, param,t)

  # a variety of g calculations depending on the model
  #g=gmax*(ii/(ii+Hi))*(N/(N+Hn))
  #g=gmax*pmin((ii/(ii+Hi)),N/(N+Hn))
   #new g
  g=gmax*pmin(Hi*ii /sqrt(gmax^2 + (Hi*ii)^2),N/(N+Hn))

  
  dPdt=-((J[2:(n+1)]-J[1:(n)])/param$dz) + (g*P)-m*P-gamma*P*P
  dNdt=-g*P+eps*Detrit-((Jdn[2:(n+1)]-Jdn[1:(n)])/param$dz) 
  # old nutrient formula.  
  #-((Jdn[2:(n+1)]-Jdn[1:(n)])/param$dz)  +alpha*m*eps*P  -alpha*g*P
  dDdt=-((JD[2:(n+1)]-JD[1:(n)])/param$dz)+m*P+gamma*P*P-eps*Detrit
  #dDdt=
  
  Y=c(dPdt,dNdt,dDdt)
  list(Y)
}

timey=1000

res=ode.1D(y=y,func=diffusivefux, times=0:timey,parms=1,nspec=1)
restid = res[,1]
resP=res[,2:(n+1)]
resN=res[,(n+2):(2*n+1)]
resD=res[,(2*n+2):(3*n+1)]




image.plot(x=restid,y=param$z,z=resP,main="concentration,depth",ylim = rev(range(param$z)),xlab="days",ylab="depth(m)")
image.plot(x=restid,y=param$z,z=resN,main="nutrients,depth",ylim = rev(range(param$z)),xlab="days",ylab="depth(m)")
image.plot(x=restid,y=param$z,z=resD,main="Detritus,depth",ylim = rev(range(param$z)),xlab="days",ylab="depth(m)")

#for ease of understanding turn off the seasonal effect 
# in the light equation
# previously to running sensitivity plots.

###### covergernce independence plot

dzrange=c(0.5,1.5,2,2.5,3)
plot(resP[(timey+1),], param$z, ylim = rev(range(param$z)), type='l', lwd=0.5, xlab = "Concentration [mmol/m^3]", ylab = "Depth [m]", main="concentration at convergence")

for (f in c(1:length(dzrange))){
  param$dz=dzrange[f]
  param$z=seq(param$dz/2,200,by=param$dz)
  n=length(param$z)
  y=c(1:(n*3))*0
  y[1]=1 # initial plankton, in nutrients.
  y[(n+1):(2*n)]=20# initial nutrients
  y[(2*n+1):(3*n)]=5# detritus.
  res=ode.1D(y=y,func=diffusivefux, times=0:timey,parms=1,nspec=1)
  restid = res[,1]
  resP=res[,2:(n+1)]
  resN=res[,(n+2):(2*n+1)]
  resD=res[,(2*n+2):(3*n+1)]
 lines( resP[(timey+1),], param$z, ylim = rev(range(param$z)), type='l',col=f+1, lwd=0.5, xlab = "Concentration [mmol/m^3]", ylab = "Depth [m]", main="concentration at convergence")
  
}
legend(0.05, 0, legend=c("dz=1", "dz=0.5", "dz=1.5", "dz=2", "dz=2.5","dz=3"),
              col=c(1:6),lwd=4, cex=0.8)
image.plot(x=restid,y=param$z,z=resP,main="concentration,depth dz=3",ylim = rev(range(param$z)),xlab="days",ylab="depth(m)")


###### sensitivity Hn

HNrange=c(0.3,3,30,300,0.03)
plot(resP[(timey+1),], param$z, ylim = rev(range(param$z)), type='l', lwd=0.5, xlab = "Concentration [mmol/m^3]", ylab = "Depth [m]", main="sensitivity analysis of Hn")

for (f in c(1:length(HNrange))){
  Hn=HNrange[f]
  res=ode.1D(y=y,func=diffusivefux, times=0:timey,parms=1,nspec=1)
  restid = res[,1]
  resP=res[,2:(n+1)]
  resN=res[,(n+2):(2*n+1)]
  resD=res[,(2*n+2):(3*n+1)]
  lines( resP[(timey+1),], param$z, ylim = rev(range(param$z)), type='l',col=f+1, lwd=4, xlab = "Concentration [mmol/m^3]", ylab = "Depth [m]", main="concentration at convergence")
  
}

legend(0.41, 0, legend=c( "Hn = 0.3", "Hn = 3","Hn = 30", "Hn = 300", "Hn = 0.03"),
       col=c(1:5),lwd=4, cex=0.8)

Hn=0.3

###### sensitivity kw
param$kw=0.1
kwrange=c(0.2,0.5,1,0.01)
plot(resP[(timey+1),], param$z, ylim = rev(range(param$z)),xlim=range(0,1.5), type='l', lwd=4, xlab = "Concentration [mmol/m^3]", ylab = "Depth [m]", main="sensitivity analysis of kw")

for (f in c(1:length(kwrange))){
  param$kw=kwrange[f]
  res=ode.1D(y=y,func=diffusivefux, times=0:timey,parms=1,nspec=1)
  restid = res[,1]
  resP=res[,2:(n+1)]
  resN=res[,(n+2):(2*n+1)]
  resD=res[,(2*n+2):(3*n+1)]
  lines( resP[(timey+1),], param$z, ylim = rev(range(param$z)), type='l',col=f+1, lwd=4, xlab = "Concentration [mmol/m^3]", ylab = "Depth [m]", main="concentration at convergence")
  
}

legend(0.8, 0, legend=c("kw=0.1","kw=0.2","kw=0.5","kw=1","kw=0.01"),
       col=c(1:5),lwd=4, cex=0.8)

