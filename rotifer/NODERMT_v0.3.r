##########################
## NODE_RMT_functions.r ##
##########################

## goal: perfrom residual minimisation training on the finch time series

## method:
# 1. interpolate the time series with sin ANN functions
# 2. estimate linear and non-linear coupling between the time series

## versions:
# 17-12-2020 - created version 0.0
# 18-12-2020 - created version 0.1
# 18-12-2020 - casted the model in a bayesian framework
# 18-12-2020 - created version 0.2
# 18-12-2020 - added a predictive model, before it was individual error models for each state variables
# 07-01-2021 - created version 0.3
# 07-01-2021 - checked code for potential errors 

## notes
# 19-12-2020 - 
# The fit is not satisfactory as the interpolation only captures the mean trend.
# The approach is supposed to simplify the NODE fitting procedure by avoiding the simulation 
# of the NODE system and by providing analytically tracktable gradients.
# Time series can be fitted individually quite easily, but that is ignoring the dynamical 
# coupling between variables, which is precisely the point of the whole endeavour.
# The method could potentially be made fast with backpropagation. 
# However, I suspect that we would still face the problem of being stuck in a local minima.
# Indeed, the interpolation of the time series greatly increases the size of the model, 
# which also decreases the capacity of the optimisation to move through the parameter space.
# I think the primary appeal for such a method was that it would potentially be much faster, and more easy to browse through the parameter space, but it turns out that it may not be the case. 
# It is certainly not suitable for fitting a "linear/polynomial" ODE... 
# From experience I know that this would be faster with a normal ODE system.
# It may be worth trying with a proper NODE system, because the possibility for computing 
# gradients may still make the approach much faster to fit with backprop.
# But even then, it would not solve the problem of the "stickiness" of the model to MCMC, 
# which is due to the larger number of parameters in the network.
# That probably makes it more suitable for fast optimisation, not really for error 
# quantification.
# To conclude, it may still be worth exploring how complicated it would be to code backprop,
# and if turns out to be to complicated, better of sticking with the NODE approach you have
# developped so far.

##############
## INITIATE ##
##############

## goal: initiate the RMT algorithm

## libraries
library(Rcpp)
sourceCpp("cpp/DEMCpp_v0.1.cpp")
sourceCpp("cpp/DEMCO_v0.5.cpp")

## TS
TS = read.table("TS.csv",sep=";",header=T)

## focus on fortis for now
TS = TS[,c("year","rain","small.seeds","large.seeds","fortis.N","fortis.PC1.beak","fortis.PC2.beak")]

## transform real variables
TS[,c("rain","small.seeds","large.seeds","fortis.N")] = log(TS[,c("rain","small.seeds","large.seeds","fortis.N")])

## test
# TS = TS[,c(1:3)]

## initiate
n     = nrow(TS)
N     = ncol(TS) - 1
sd1_o = 0.1
sd2_o = 0.5
sd1_p = 0.1
sd2_p = 0.5
alpha = 10
t     = TS[,1]
Y     = TS[,-1] 
t     = (t-min(t))/(max(t)-min(t))*alpha
mean_ = apply(Y,2,mean)
sd_   = apply(Y,2,sd)
Y     = t(apply(Y,1,function(x)(x-mean_)/sd_)) # standardise Y
t_    = seq(min(t),max(t),(max(t)-min(t))/100)

## functions for the observation model
f_o = function(t,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	return(t(Omega[,1])%*%sin(pi*(t*Omega[,2] + Omega[,3])))
}
ddt.f_o = function(t,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	return(t(Omega[,1])%*%(pi*Omega[,2]*cos(pi*(t*Omega[,2] + Omega[,3]))))
}
f_o.eval     = function(t,Omega) apply(t(t),2,function(x)f_o(x,Omega))
ddt.f_o.eval = function(t,Omega) apply(t(t),2,function(x)ddt.f_o(x,Omega))
ddOmega_o.f_o = function(x,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	dfdOmega1 = sin(pi*(x*Omega[,2] + Omega[,3]))
	dfdOmega2 = Omega[,1]*pi*x*cos(pi*(x*Omega[,2] + Omega[,3]))
	dfdOmega3 = Omega[,1]*pi*1*cos(pi*(x*Omega[,2] + Omega[,3]))
	return(c(dfdOmega1,dfdOmega2,dfdOmega3))
}
ddOmega_o.f_o.eval = function(x,Omega)
{
	apply(t(x),2,function(x) ddOmega_o.f_o(x,Omega))
}
ddOmega_o.ddt.f_o = function(x,Omega) # not safe for work
{	
	Omega = matrix(Omega,ncol=3)
	dfdOmega1 = pi*Omega[,2]*cos(pi*(x*Omega[,2] + Omega[,3]))
	dfdOmega2 = Omega[,1]*(pi*cos(pi*(x*Omega[,2] + Omega[,3])) - pi*Omega[,2]*x*sin(pi*(x*Omega[,2] + Omega[,3])))
	dfdOmega3 = -Omega[,1]*(pi*Omega[,2]*1*sin(pi*(x*Omega[,2] + Omega[,3])))
	return(c(dfdOmega1,dfdOmega2,dfdOmega3))
}

## functions for the process model
# PIM = function(x) # pairwise interaction matrix function
# { 
# 	PIM_         = x%*%t(x)
# 	PIM_diag     = diag(PIM_)
# 	PIM_lowertri = as.vector(PIM_)[as.vector(lower.tri(PIM_))]
# 	return(c(PIM_diag,PIM_lowertri))
# }
ddt.f_p           = function(x,Omega) x*Omega
ddt.f_p.eval      = function(X,Omega) X%*%Omega
ddOmega_p.ddt.f_p = function(X,Omega) X

predict = function(t,t_,Omega_o,Omega_p) 
{
	## initiate
	Omega_o = matrix(Omega_o,ncol=N)
	Omega_p = matrix(Omega_p,ncol=N)


	## observation model 
	Ybar_o      = apply(Omega_o,2,function(x) f_o.eval(t,x))
	Ybar_o_     = apply(Omega_o,2,function(x) f_o.eval(t_,x))
	ddt.Ybar_o_ = apply(Omega_o,2,function(x) ddt.f_o.eval(t_,x))

	## process model
	X_p = cbind(1,Ybar_o_)
	ddt.Ybar_p_ = apply(Omega_p,2,function(x) ddt.f_p.eval(X_p,x))

	return(list("Ybar_o"=Ybar_o,"Ybar_o_"=Ybar_o_,"ddt.Ybar_o_"=ddt.Ybar_o_,"ddt.Ybar_p_"=ddt.Ybar_p_))
}

logPost = function(t,t_,Omega_o,Omega_p)
{
	## initiate
	Omega_o = matrix(Omega_o,ncol=N)
	Omega_p = matrix(Omega_p,ncol=N)

	## observation model 
	Ybar_o      = apply(Omega_o,2,function(x) f_o.eval(t,x))
	Ybar_o_     = apply(Omega_o,2,function(x) f_o.eval(t_,x))
	ddt.Ybar_o_ = apply(Omega_o,2,function(x) ddt.f_o.eval(t_,x))

	## process model
	X_p = cbind(1,Ybar_o_)
	ddt.Ybar_p_ =  apply(Omega_p,2,function(x) ddt.f_p.eval(X_p,x))
	
	## residuals
	res_o = (Y - Ybar_o)
	res_p = (ddt.Ybar_o_ - ddt.Ybar_p_)

	## compute posterior
	logLik_o   = sum(log(dnorm(res_o,0,sd1_o)))
	logLik_p   = sum(log(dnorm(res_p,0,sd1_p)))
	logPrior_o = sum(log(dnorm(Omega_o,0,sd2_o)))
	logPrior_p = sum(log(dnorm(Omega_p,0,sd2_p)))
	logPost    = 1.0*(logLik_o + logPrior_o) + 0.0*(logLik_p + logPrior_p)

	## terminate
	return(logPost)
}

logPost2 = function(t,t_,Omega_o,Omega_p)
{
	## initiate
	Omega_o = matrix(Omega_o,ncol=N)
	Omega_p = matrix(Omega_p,ncol=N)

	## observation model 
	Ybar_o      = apply(Omega_o,2,function(x) f_o.eval(t,x))
	Ybar_o_     = apply(Omega_o,2,function(x) f_o.eval(t_,x))
	ddt.Ybar_o_ = apply(Omega_o,2,function(x) ddt.f_o.eval(t_,x))

	## process model
	X_p = cbind(1,Ybar_o_)
	ddt.Ybar_p_ =  apply(Omega_p,2,function(x) ddt.f_p.eval(X_p,x))
	
	## residuals
	res_o = (Y - Ybar_o)
	res_p = (ddt.Ybar_o_ - ddt.Ybar_p_)

	## compute posterior
	logLik_o   = sum(-0.5*log(2*pi*sd1_o^2) - 0.5*(res_o^2)/sd1_o^2)
	logLik_p   = sum(-0.5*log(2*pi*sd1_p^2) - 0.5*(res_p^2)/sd1_p^2)
	logPrior_o = sum(-0.5*log(2*pi*sd2_o^2) - 0.5*(Omega_o^2)/sd2_o^2)
	logPrior_p = sum(-0.5*log(2*pi*sd2_p^2) - 0.5*(Omega_p^2)/sd2_p^2)
	logPost    = logLik_o + logLik_p + logPrior_o + logPrior_p

	## terminate
	return(logPost)
}

ddOmega.logPost2 = function(t,t_,Omega_o,Omega_p)
{
	## initiate
	Omega_o = matrix(Omega_o,ncol=N)
	Omega_p = matrix(Omega_p,ncol=N)

	## observation model 
	Ybar_o      = apply(Omega_o,2,function(x) f_o.eval(t,x))
	Ybar_o_     = apply(Omega_o,2,function(x) f_o.eval(t_,x))
	ddt.Ybar_o_ = apply(Omega_o,2,function(x) ddt.f_o.eval(t_,x))

	## process model
	X_p = cbind(1,Ybar_o_)
	ddt.Ybar_p_ =  apply(Omega_p,2,function(x) ddt.f_p.eval(X_p,x))
	
	## residuals
	res_o = Y - Ybar_o
	res_p = ddt.Ybar_o_ - ddt.Ybar_p_
	ddOmega.res_o = -ddOmega.Ybar_o
	ddOmega.res_p = ddOmega.ddt.Ybar_o_ - ddOmega.ddt.Ybar_p_

	## compute posterior
	ddOmega.logLik_o   = -sum(ddOmega.res_o*res_o/sd1_o^2)
	ddOmega.logLik_p   = -sum(ddOmega.res_p*res_p/sd1_p^2)
	ddOmega.logPrior_o = -sum(Omega_o/sd2_o^2)
	ddOmega.logPrior_p = -sum(Omega_p/sd2_p^2)
	ddOmega.logPost    = ddOmega.logLik_o + ddOmega.logLik_p + ddOmega.logPrior_o + ddOmega.logPrior_p

	## terminate
	return(logPost)
}

#
###

###############
## FIT MODEL ##
###############

## goal: fit the model to the time series

## optimise model
## BFGS
# Omega  = rnorm(N*3*10+N*(N+1),0,0.1)
# error_ = function(x)-logPost(t,t_,x[1:(N*3*10)],x[-c(1:(N*3*10))])
# Omega  = optim(par=Omega,fn=error_,method="BFGS",control=list("trace"=1,"REPORT"=1,"maxit"=10))$par

## DEMC
# error_ = function(x)logPost(t,t_,x[1:(N*3*10)],x[-c(1:(N*3*10))])
# Omega  = rnorm(N*3*10+N*(N+1),0,0.1)
# chainList = DEMCOpp(list("dTarget"=error_,"Theta_0"=Omega,"epsilon"=0.001,"lambda"=10,"nIt"=1000))
# Omega = chainList$chain[which.max(chainList$chain[,1]),-1]

## DEMC
# error_ = function(x)logPost(t,t_,x[1:(N*3*10)],x[-c(1:(N*3*10))])
# Omega  = rnorm(N*3*10+N*(N+1),0,0.1)
# chainList = DEMCpp(list("dTarget"=error_,"Theta_0"=Omega,"epsilon"=0.001,"nIt"=1000))
# Omega = chainList$chain[which.max(chainList$chain[,1]),-1]

#
###

##############
## ANALYSIS ##
##############

pdf("results.pdf")

## compute predictions
Omega_o = Omega[1:(N*3*10)]
Omega_p = Omega[-c(1:(N*3*10))]
predicted = predict(t,t_,Omega_o,Omega_p)
Ybar_o  = predicted$Ybar_o_
dYbar_o = predicted$ddt.Ybar_o_
dYbar_p = predicted$ddt.Ybar_p_

## visualise fit
par(mfrow=c(3,1))
colVect = rainbow(N)
for(i in 1:N)
{
	plot(t,Y[,i],type="b",ylim=c(-3,3))
	lines(t_,Ybar_o[,i],col=rainbow(N)[i])
	lines(c(min(t),max(t)),c(0,0),lty=2)
	#
	plot(t_,dYbar_o[,i],col="blue",ylim=c(-4,4),type="l")
	lines(t_,dYbar_p[,i],col="green")
	lines(c(min(t),max(t)),c(0,0),lty=2)
	#
	Omega_p = matrix(Omega_p,ncol=N)
	plot(t_,rep(Omega_p[1,i],length(t_)),ylim=c(-max(Omega_p[,i]),max(Omega_p[,i])),type="l")
	for(j in 1:N){lines(t_,rep(Omega_p[1+j,i],length(t_)),col=rainbow(N)[j])}
	lines(c(min(t),max(t)),c(0,0),lty=2)
}
par(mfrow=c(1,1))

## visualise phase space

# ## residuals checks
# s = match(x=round(t,2),table=round(t_,2))
# par(mfrow=c(2,2))
# for(i in 1:N)
# {
# 	res_o = Y - Ybar_o[s,]
# 	res_p = dYbar_o - dYbar_p 
# 	barplot(res_o[,i])
# 	barplot(res_p[,i])
# }
# par(mfrow=c(1,1))

dev.off()

#
###
