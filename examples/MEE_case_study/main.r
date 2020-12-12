############
## main.r ##
############

## goal: code supporting the case study of the 2020 MEE paper "Neural ordinary differential equations for 
# ecological and evolutionary time series anlaysis" W. Bonnaffe, B.C. Sheldon, T. Coulson

## update
## 12-12-2020 - created file

##############
## INITIATE ##
##############

## import NODE functions
source("https://raw.githubusercontent.com/WillemBonnaffe/NODER/master/versions/v1.0/NODER_v1.0.0.r")

#
###

###################
## NODE INITIATE ##
###################

## load time series
TS = read.table("TS.csv",sep=";",header=T)

## format time series
colnames(TS) = c("t","H","L") # compact names
TS_ = TS.format(TS,c(2,3))    # log transform and standardise

## initiate model
M0 = NODE.init(TS=TS_,
			   SLPinputStateList=list(c(2,3),c(2,3)),
			   SLPftypeVect=c(2,2)
			   )

## or load previous model
# load("model_P0.15.RData")
# load("model_P0.3.RData")

## user defined parameters for the normal Bayesian model
sd_lik = rep(0.1,M0$N-1)
sd_prior_network = rep(0.3,M0$N-1)
#
## initiate sd_lik
M0$sd_lik = rep(sd_lik,rep(M0$n,M0$N-1))
#
## initiate mu prior
mu_prior_states  = as.vector(t(M0$TS[1,-1])) # set mu_prior_states to first observation
mu_prior_network = rep(0,M0$D-(M0$N-1))
M0$mu_prior         = c(mu_prior_states,mu_prior_network)
#
## initiate sd prior
sd_prior_states  = as.vector(t(apply(M0$TS[,-1],2,sd))) # set sd_prior_states to sd of time series
sd_prior_network = rep(sd_prior_network,M0$dimVect[-1]) 
M0$sd_prior         = c(sd_prior_states,sd_prior_network)

#
###

##############
## NODE FIT ##
##############

# ## prefit with LASSO/Ridge regression
# # of dYdt operator on raw time differences from time series
# error.prefit = function(Omega)
# {
#   .error.LR(Y      =  TS.timedifferentiate(M0$TS)[,-1],
#             X      =  M0$TS[-M0$n,],
#             f      =  function(X,param){t(apply(X,1,function(x)M0$.Ydot(x,param)[-1]))},
#             param  = Omega,
#             lambda = 0.002,
#             alpha  = 1
#             )[[1]] 
# }
# Omega = rnorm(M0$D-(M0$N-1),0,0.02); print(error.prefit(Omega))
# optList = .fit(Omega,error.prefit,20)
# Omega = optList$par

## prefit with normal bayesian
# of dYdt operator on raw time differences from time series
error.prefit = function(Omega)
{
  -.error.normalBayesian(Y        = TS.timedifferentiate(M0$TS)[,-1],
                        X         = M0$TS[-M0$n,],
                        f         = function(X,param){t(apply(X,1,function(x)M0$.Ydot(x,param)[-1]))},
                        param     = Omega,
                        sd_lik    = M0$sd_lik,
                        mu_prior  = M0$mu_prior,
                        sd_prior  = M0$sd_prior
                        )[[1]] 
}
Omega = rnorm(M0$D-(M0$N-1),0,0.001); print(error.prefit(Omega))
optList = .fit(Omega,error.prefit,25)
Omega = optList$par

## fit
error.fit = function(Theta)
{
  -.error.normalBayesian(Y         = M0$TS[,-1],
                         X         = M0$TS[,1],
                         f         = function(X,param){M0$.Ybar.ode(X,param[M0$i_init],param[-M0$i_init])[,-1]},
                         param     = Theta,
                         sd_lik    = M0$sd_lik,
                         mu_prior  = M0$mu_prior,
                         sd_prior  = M0$sd_prior
                         )[[1]]
}
Theta = c(rep(0,M0$N-1),Omega) # starting from pre-fit
# Theta = rnorm(M0$D,0,0.01) # random starting point
print(error.fit(Theta)) # check
optList = .fit(Theta,error.fit,100)
M0$Theta = Theta = optList$par

## estimate error with ABC
chainList = .ABC(M0$Theta,error.fit,1000,noise=0.01,threshold=0.95)
M0$chain = chainList$chain

## save model
save(M0,file="model.RData")

#
###

###################
## NODE ANALYSIS ##
###################

pdf("results.pdf")
par(cex.lab=1.5,oma=c(1,1,1,1))

## visualise dynamics
t = seq(min(M0$TS[,1]),max(M0$TS[,1]),.1)
Y_0 = M0$Theta[M0$i_init]
Omega = M0$Theta[-M0$i_init]
Ybar = M0$.Ybar.ode(t,Y_0,Omega)
.plot.fit(M0$TS,Ybar)

## compute approximated posterior distribution
t = seq(0,90,1)
#
## states
Ybar.ensemble = apply(M0$chain[,-1],1,function(x)M0$.Ybar.ode(t,x[M0$i_init],x[-M0$i_init]))
Ybar.q05 = matrix(apply(Ybar.ensemble,1,function(x)quantile(x,probs=c(0.05))),ncol=M0$N)
Ybar.q50 = matrix(apply(Ybar.ensemble,1,function(x)quantile(x,probs=c(0.50))),ncol=M0$N)
Ybar.q95 = matrix(apply(Ybar.ensemble,1,function(x)quantile(x,probs=c(0.95))),ncol=M0$N)
#
## jacobian
Jbar.ensemble = apply(M0$chain[,-1],1,function(x)M0$.Jbar.ode(t,x[M0$i_init],x[-M0$i_init]))
Jbar.q05 = matrix(t(apply(Jbar.ensemble,1,function(x)quantile(x,probs=c(0.05)))),nrow=length(t))
Jbar.q50 = matrix(t(apply(Jbar.ensemble,1,function(x)quantile(x,probs=c(0.50)))),nrow=length(t))
Jbar.q95 = matrix(t(apply(Jbar.ensemble,1,function(x)quantile(x,probs=c(0.95)))),nrow=length(t))
#
## contributions
Cbar.ensemble = apply(M0$chain[,-1],1,function(x)M0$.Cbar.ode(t,x[M0$i_init],x[-M0$i_init]))
Cbar.q05 = matrix(t(apply(Cbar.ensemble,1,function(x)quantile(x,probs=c(0.05)))),nrow=length(t))
Cbar.q50 = matrix(t(apply(Cbar.ensemble,1,function(x)quantile(x,probs=c(0.50)))),nrow=length(t))
Cbar.q95 = matrix(t(apply(Cbar.ensemble,1,function(x)quantile(x,probs=c(0.95)))),nrow=length(t))

## plot fit, effects, contributions
par(mar=c(4,4.1,2,2))
layout(mat=rbind(c(1,4),c(2,5),c(3,6)))
colVect = rainbow(M0$N)
#
for(i in 2:M0$N)
{
	## states
	X.q05 = Ybar.q05
	X.q50 = Ybar.q50
	X.q95 = Ybar.q95
	plot(M0$TS[,1],M0$TS[,i],col=colVect[i],xlim=c(min(t),max(t)),ylim=c(-1,1)*3,xlab="t",ylab=colnames(M0$TS)[i],pch=16)
	polygon(x=c(t,rev(t)),y=c(X.q05[,i],rev(X.q95[,i])),border=NA,col=adjustcolor(colVect[i],alpha=0.2))
	lines(t,X.q50[,i],col=colVect[i])
	#
	## jacobian
	X.q05 = Jbar.q05
	X.q50 = Jbar.q50
	X.q95 = Jbar.q95
	plot(x=c(min(t),max(t)),y=c(0,0),lty=2,ylim=c(-1,1)*1.5,cex=0,xlab="t",ylab="Effects")
	legend("bottom",col=colVect,legend=colnames(M0$TS),lty=1,bty="n",horiz=T,cex=1.5)
	for(j in 1:M0$N)
	{
		k = (j-1)*(M0$N) + i
		polygon(x=c(t,rev(t)),y=c(X.q05[,k],rev(X.q95[,k])),border=NA,col=adjustcolor(colVect[j],alpha=0.2))
		lines(t,X.q50[,k],col=colVect[j])
	}
	#
	## contributions
	X.q05 = Cbar.q05
	X.q50 = Cbar.q50
	X.q95 = Cbar.q95
	plot(x=c(min(t),max(t)),y=c(0,0),lty=2,ylim=c(-1,1)*1.0,cex=0,xlab="t",ylab="Contributions")
	for(j in 1:M0$N)
	{
		k = (j-1)*(M0$N) + i
		polygon(x=c(t,rev(t)),y=c(X.q05[,k],rev(X.q95[,k])),border=NA,col=adjustcolor(colVect[j],alpha=0.2))
		lines(t,X.q50[,k],col=colVect[j])
	}
}
#
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1)

par(mfrow=c(2,2))
#
## visualise dYdt functions
for(i in c(2,3))
{
	dYdt = function(x){M0$.Ydot(c(0,x[1],x[2]),M0$Theta[-M0$i_init])[i]}
	.plotMat(lims=rbind(c(-2,2),c(-2,2)),func=dYdt)
	TS_ = (M0$TS-(-2))/(2-(-2))
	points(TS_[,2],TS_[,3],pch=16)
}
#
## average Jacobian and contributions
Jbar.q50.mean = matrix(apply(Jbar.q50,2,mean),ncol=M0$N)[-1,-1]
Cbar.q50.mean = matrix(apply(abs(Cbar.q50),2,mean),ncol=M0$N)[-1,-1]
#
## dynamical interaction network (DIN) plot
.plot.DIN(Jbar.q50.mean,Cbar.q50.mean,labels=colnames(M0$TS)[-1])
#
## pairwise phase space (PPS) plot
dYdt = function(x){M0$.Ydot(c(0,x[1],x[2]),M0$Theta[-M0$i_init])}
.plot.PPS(2,3,c(2,3),dYdt,rbind(c(-2,2),c(-2,2)))
#
## find equilibrium
Y_eq = rnorm(M0$N-1,0,1)
dYdt = function(x){sum((M0$.Ydot(c(0,x),M0$Theta[-M0$i_init]))^2)}
Y_eq = .fit(Y_eq,dYdt,100)
Y_eq_scaled = (Y_eq$par-(-2))/(2-(-2))
points(Y_eq_scaled[1],Y_eq_scaled[2],col="red",pch=8,cex=2)
#
par(mfrow=c(1,1))

## check residuals
res = M0$TS[,-1] - M0$.Ybar.ode(M0$TS[,1],Y_0,Omega)[,-1]
.plot.residualCheck1(res)
.plot.residualCheck2(res)
.plot.residualCheck3(res)
r2 = 1-sd(t(res))^2/sd(t(M0$TS[,-1]))^2

par(cex.lab=1)
dev.off()

#
###
