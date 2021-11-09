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

##############
## INITIATE ##
##############

## goal: initiate the RMT algorithm

## TS
TS = read.table("TS.csv",sep=";",header=T)

## focus on fortis for now
TS = TS[,c("year","rain","small.seeds","large.seeds","fortis.N","fortis.PC1.beak","fortis.PC2.beak")]

## transform real variables
TS[,c("rain","small.seeds","large.seeds","fortis.N")] = log(TS[,c("rain","small.seeds","large.seeds","fortis.N")])

## initiate
n     = nrow(TS)
N     = ncol(TS) - 1
sd1_o = 0.5
sd2_o = 0.5
sd1_p = 0.5
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
f_o.eval     = function(t,Omega) apply(t(t),2,function(x)f_o(t,Omega))
ddx.f_o.eval = function(t,Omega) apply(t(t),2,function(x)ddt.f_o(t,Omega))
error.f_o = function(Y,Omega)
{
	res      = Y-f_o.eval(t,Omega)
	logLik   = sum(log(dnorm(res,0,sd1_o)))
	logPrior = sum(log(dnorm(Omega,0,sd2_o)))
	logPost  = logLik + logPrior
	return(logPost)
}
# ddOmega_o.f_o = function(x,Omega)
# {	
# 	Omega = matrix(Omega,ncol=3)
# 	dfdOmega1 = sin(pi*(x*Omega[,2] + Omega[,3]))
# 	dfdOmega2 = Omega[,1]*pi*x*cos(pi*(x*Omega[,2] + Omega[,3]))
# 	dfdOmega3 = Omega[,1]*pi*1*cos(pi*(x*Omega[,2] + Omega[,3]))
# 	return(c(dfdOmega1,dfdOmega2,dfdOmega3))
# }
# ddOmega_o.f_o.eval = function(x,Omega)
# {
# 	apply(t(x),2,function(x)ddOmega_o.f_o(x,Omega))
# }
# ddOmega_o.error.f_o = function(Y,Omega){as.vector(1/n*2*(-ddOmega_o.f_o.eval(t,Omega))%*%(Y-f_o.eval(t,Omega)) + lambda_p*2*Omega)} # not safe for work

## functions for the process model
# PIM = function(x) # pairwise interaction matrix function
# { 
# 	PIM_         = x%*%t(x)
# 	PIM_diag     = diag(PIM_)
# 	PIM_lowertri = as.vector(PIM_)[as.vector(lower.tri(PIM_))]
# 	return(c(PIM_diag,PIM_lowertri))
# }
X_p          = cbind(1,Ybar_o) #,t(apply(Ybar_o,1,PIM)))
ddt.f_p      = function(x,Omega) x*Omega
ddt.f_p.eval = function(Ybar_o,Omega)
{	
	apply(Ybar_o,1,Omega)
}
error.f_p    = function(ddt.f_o,Omega) 
{
	res      = ddx.f_o - ddx.f_p(X_p,Omega)
	logLik   = sum(log(dnorm(res,0,sd1_p)))
	logPrior = sum(log(dnorm(Omega,0,sd2_p)))
	logPost  = logLik + logPrior
	return(logPost)
}
# ddOmega_p.ddx.f_p   = function(X,Omega) X
# ddOmega_p.error.f_p = function(ddx.f_o,Omega) 1/n*2*t(-ddOmega_p.ddx.f_p(X_p,Omega))%*%(ddx.f_o - ddx.f_p(X_p,Omega)) + lambda_p*2*Omega # not safe for work

predict = function(t,t_,Omega_o,Omega_p) 
{
	## initiate
	Omega_o = matrix(Omega_o,ncol=N)
	Omega_p    = matrix(Omega_p,ncol=N)

	## evaluate
	Ybar_o      = apply(Omega_o,2,function(x) f_o.eval(t,x))
	Ybar_o_     = apply(Omega_o,2,function(x) f_o.eval(t_,x))
	ddt.Ybar_o_ = apply(Omega_o,2,function(x) ddt.f_o.eval(t_,x))
	ddt.Ybar_p_ = apply(Omega_p,2,function(y) apply(Ybar_o_,1,ddx.f_p(x,y)))
}

## posterior distribution
error = function(Omega)
{
	error_  = 0
	Omega_o = matrix(Omega[1:(N*3*10)],ncol=N)
	Omega_p = matrix(Omega[-c(1:(N*3*10))],ncol=N)
	for(i in 1:N)
	{
		error_ = error_ + 
	 		 	 error.f_o(Y[,i],Omega_o[,i]) + 
			 	 error.f_p(ddx.f_o.eval(t_,Omega_o[,i]),Omega_p[,i])
	}
	return(error_)
}

#
###

###############
## FIT MODEL ##
###############

## goal: fit the model to the time series

# ## optimise model
# Omega = rnorm(N*3*10+N*ncol(X_p),0,0.1)
Omega = optim(par=Omega,fn=function(x)-error(x),method="BFGS",control=list("trace"=1,"REPORT"=1,"maxit"=10))$par

#
###

##############
## ANALYSIS ##
##############

pdf("results.pdf")

## compute predictions
Omega_o = matrix(Omega[1:(N*3*10)],ncol=N)
Omega_p = matrix(Omega[-c(1:(N*3*10))],ncol=N)

## fit observation model
Ybar_o  = NULL
dYbar_o = NULL
for(i in 1:N)
{
	## predict
	Ybar_o_  = f_o.eval(t_,Omega_o[,i])
	dYbar_o_ = ddx.f_o.eval(t_,Omega_o[,i])

	## store predictions
	Ybar_o  = cbind(Ybar_o,Ybar_o_)
	dYbar_o = cbind(dYbar_o,dYbar_o_)
}
# plot(data.frame(Ybar_o))

## fit process model
dYbar_p = NULL
for(i in 1:N) 
{
	## predict
	dYbar_p_    = ddx.f_p(X_p,Omega_p[,i])

	## store predicition
	dYbar_p     = cbind(dYbar_p,dYbar_p_)
}


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
	# lines(t_,dYbar_p_num[,i],col="orange") # for debugging
	#
	plot(t_,rep(Omega_p[1,i],length(t_)),ylim=c(-max(Omega_p[,i]),max(Omega_p[,i])),type="l")
	for(j in 1:N){lines(t_,rep(Omega_p[1+j,i],length(t_)),col=rainbow(N)[j])}
	lines(c(min(t),max(t)),c(0,0),lty=2)
}
par(mfrow=c(1,1))

## visualise phase space

## residuals checks
s = match(x=round(t,2),table=round(t_,2))
par(mfrow=c(2,2))
for(i in 1:N)
{
	res_o = Y - Ybar_o[s,]
	res_p = dYbar_o - dYbar_p 
	barplot(res_o[,i])
	barplot(res_p[,i])
}
par(mfrow=c(1,1))

dev.off()

#
###
