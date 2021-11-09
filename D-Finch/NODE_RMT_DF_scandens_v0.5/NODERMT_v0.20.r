##########################
## NODE_RMT_functions.r ##
##########################

## goal: 
# - perfrom residual minimisation training on the time series data
# - analyse the temporal change in variables with NODEs

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
# 07-01-2021 - created version 0.4
# 07-01-2021 - replaced polynomial process model by ANN process model 
# 08-01-2021 - created version 0.5
# 08-01-2021 - decoupled the interpolation from the fitting of the process model
# 14-01-2021 - implemented visualisation of the interpolation 
# 15-01-2021 - created version 0.6
# 15-01-2021 - implemented a bifurcation on the level of complexity considered
# 15-01-2021 - removed DEMCO as now obsolete because using gradients
# 15-01-2021 - computed r2 for each time series
# 15-01-2021 - created version 0.7
# 15-01-2021 - implemented computation of identifiability by ensembling
# 25-01-2021 - implemented computation of BIC of each model to formalise identification of complexity level 
# 30-01-2021 - implemented computation of the minimum acceptable prior complexity
# 03-02-2021 - created version 0.8
# 03-02-2021 - made sure that code was applicable to other time series by applying it to the tri-trophic system (Algee, Flagellate, Rotifer)
# 05-02-2021 - created version 0.9
# 05-02-2021 - implemented the bayesian regularisation scheme developped by Cawley and Talbot
# 14-02-2021 - created version 0.10
# 14-02-2021 - debugged code, switched back to bayesian regularisation with prior variance
# 15-02-2021 - implemented selection on complexity by minimising variance around estimates of effects through visual assessments
# 15-02-2021 - simplified code
# 22-02-2021 - created version 0.11
# 22-02-2021 - RE-implemented the Cawley approach
# 22-02-2021 - found a stable combination of regularisation parameters to fit the finch series with the standard normal approach
# 23-02-2021 - created version 0.12
# 23-02-2021 - tried to make the Cawley approach work 
# 24-02-2021 - created version 0.13
# 24-02-2021 - tried again the maximum constraint approach
# 24-02-2021 - created version 0.14
# 24-02-2021 - tested MCMC sampling instead to find the expected interpolation 
#            => it is not better, the model is still not fitted after 1,000,000 iteration
# 24-02-2021 - created version 0.15
#            - removed MCMC sampling
# 24-02-2021 - created version 0.16
#            - settled for the standard normal approach to regularisation (results will require a sensitivity analysis on prior sd)
#            - compared results coming from raw differences
# 24-02-2021 - created version 0.16
#            - settled for the standard normal approach to regularisation (results will require a sensitivity analysis on prior sd)
#            - compared results coming from raw differences
# 05-03-2021 - created version 0.17
#            - applied to the rotifer time series again
# 09-03-2021 - created version 0.18
#            - cleaned the code
#            - made code applicable to any time series
# 15-03-2021 - created version 0.19
#            - made figures for paper
# 17-03-2021 - polished figures
# 29-03-2021 - created version 0.20
#            - separated the fitting and the analysis/plotting

#################################
## FUNCTIONS OBSERVATION MODEL ##
#################################

## functions for the observation model

## f_o
## goal: interpolated state variable
f_o = function(t,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	return(t(Omega[,1])%*%sin(pi*(t*Omega[,2] + Omega[,3])))
}

## ddt.f_o
## goal: time derivative of the interpolated variable
ddt.f_o = function(t,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	return(t(Omega[,1])%*%(pi*Omega[,2]*cos(pi*(t*Omega[,2] + Omega[,3]))))
}

## ddOmega_o.f_o
## goal: derivative of the interpolation wtr to each network parameter
ddOmega_o.f_o = function(x,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	dfdOmega1 = sin(pi*(x*Omega[,2] + Omega[,3]))
	dfdOmega2 = Omega[,1]*pi*x*cos(pi*(x*Omega[,2] + Omega[,3]))
	dfdOmega3 = Omega[,1]*pi*1*cos(pi*(x*Omega[,2] + Omega[,3]))
	return(c(dfdOmega1,dfdOmega2,dfdOmega3))
}

## ddOmega_o.ddt.f_o
## goal: derivative of the time derivative of the interpolation wtr to each network parameter
ddOmega_o.ddt.f_o = function(x,Omega) 
{	
	Omega = matrix(Omega,nrow=W_o)
	ddOmega_o1. = pi*Omega[,2]*cos(pi*(x*Omega[,2] + Omega[,3]))
	ddOmega_o2. = Omega[,1]*(pi*cos(pi*(x*Omega[,2] + Omega[,3])) - pi*Omega[,2]*x*sin(pi*(x*Omega[,2] + Omega[,3])))
	ddOmega_o3. = -Omega[,1]*(pi*Omega[,2]*1*sin(pi*(x*Omega[,2] + Omega[,3])))
	return(c(ddOmega_o1.,ddOmega_o2.,ddOmega_o3.))
}

## evaluate functions across multiple time steps
## goal: evaluate functions across multiple times
f_o.eval           = function(t,Omega) apply(t(t),2,function(x) f_o(x,Omega))
ddt.f_o.eval       = function(t,Omega) apply(t(t),2,function(x) ddt.f_o(x,Omega))
ddOmega_o.f_o.eval = function(t,Omega) apply(t(t),2,function(x) ddOmega_o.f_o(x,Omega))

## optimisation on the first level of inference
# #
# logPost_o = function(t,Y,Omega_o)
# {
# 	## predict
# 	Ybar_o = f_o.eval(t,Omega_o)
# 
# 	## residuals
# 	res_o = Y - Ybar_o
# 
# 	## Likelihood
# 	logLik_o = - sum(res_o^2)/(sd1_o^2)
# 
# 	## prior
# 	logPrior_o = - sum(Omega_o^2)/(sd2_o^2)
# 
# 	## posterior
# 	logPost_o  = logLik_o + logPrior_o
# 
# 	## terminate
# 	return(logPost_o)
# }
# # 
# ddOmega_o.logPost_o = function(t,Y,Omega_o)
# {
# 	## predict
# 	Ybar_o = f_o.eval(t,Omega_o)
# 
# 	## derivative of residuals
# 	res_o = (Y - Ybar_o)
# 	ddOmega_o.res_o = - ddOmega_o.f_o.eval(t,Omega_o)
# 
# 	## derivative of likelihood
# 	ddOmega_o.logLik_o = - 2 * ddOmega_o.res_o%*%res_o/(sd1_o^2)
# 
# 	## derivative of prior
# 	ddOmega_o.logPrior_o = - 2 * Omega_o/(sd2_o^2)
# 
# 	## derivative of posterior
# 	ddOmega_o.logPost_o = ddOmega_o.logLik_o + ddOmega_o.logPrior_o
# 
# 	## terminate
# 	return(ddOmega_o.logPost_o)
# }

## optimisation on the second level of inference
#
logPost_o = function(t,Y,Omega_o)
{
	## predict
	Ybar_o = f_o.eval(t,Omega_o)

	## residuals
	res_o = Y - Ybar_o

	## Likelihood
	logLik_o = - 0.5 * length(Y) * log(0.5 * sum(res_o^2) + 1)

	## prior
	logPrior_o = - 0.5 * length(Omega_o) * log(0.5 * sum(Omega_o^2) + 1)

	## posterior
	logPost_o  = logLik_o + logPrior_o

	## terminate
	return(logPost_o)
}
#
ddOmega_o.logPost_o = function(t,Y,Omega_o)
{
	## predict
	Ybar_o = f_o.eval(t,Omega_o)

	## derivative of residuals
	res_o = (Y - Ybar_o)
	ddOmega_o.res_o = -ddOmega_o.f_o.eval(t,Omega_o)

	## derivative of likelihood
	ddOmega_o.logLik_o = - 0.5 * length(Y) * 1/(0.5 * sum(res_o^2) + 1) * 0.5 * ddOmega_o.res_o%*%res_o

	## derivative of prior
	ddOmega_o.logPrior_o = - 0.5 * length(Omega_o) * 1/(0.5 * sum(Omega_o^2) + 1) * Omega_o

	## derivative of posterior
	ddOmega_o.logPost_o = ddOmega_o.logLik_o + ddOmega_o.logPrior_o

	## terminate
	return(ddOmega_o.logPost_o)
}

## argmax.logPost_o
# goal: maximise the log posterior density of the parameters
argmax.logPost_o = function(t,Y,Omega_o)
{
	error_     = function(x) -logPost_o(t,Y,x)
	graderror_ = function(x) -ddOmega_o.logPost_o(t,Y,x)
	Omega_o  = optim(par=Omega_o,
					 fn=error_,
					 gr=graderror_,
					 method="BFGS"# ,
					 # control=list("trace"=1,"REPORT"=1,"maxit"=100)
					 )$par
	return(Omega_o)
}

#
###

#############################
## FUNCTIONS PROCESS MODEL ##
#############################

## goal: define functions to fit the process model to the time series

## activation functions
f_sigma     = function(x) 1/(1+exp(-x))
ddu.f_sigma = function(x) f_sigma(x) * (1 - f_sigma(x))
# f_sigma     = function(x) x
# ddu.f_sigma = function(x) 1

## ddt.f_p 
# goal: process model defined as a neural ordinary differential equation (NODE) 
#       to approximate dynamics of each state variable
ddt.f_p = function(x,Omega) 
{	
	Omega = matrix(Omega,nrow=W_p)
	return(t(Omega[,1])%*%f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x))
}

## ddOmega_p.ddt.f_p
# goal: derivative of the process model wtr to parameters
ddOmega_p.ddt.f_p = function(x,Omega)
{
	Omega = matrix(Omega,nrow=W_p)
	ddOmega_p1. = f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x)
	ddOmega_p2. = Omega[,1]*1*ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x)
	ddOmega_p3. = Omega[,1]*rep(1,W_p)%*%t(x)*as.vector(ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x))
	return(c(ddOmega_p1.,ddOmega_p2.,ddOmega_p3.))
}
#
## ddx.ddt.f_p
# goal: derivative of the process model wtr to state variables
ddx.ddt.f_p = function(x,Omega)
{
	Omega = matrix(Omega,nrow=W_p)
	ddx. = Omega[,1]%*%(Omega[,-c(1:2)]*as.vector(ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x)))
	return(ddx.)
}

## evaluate process functions across several timesteps
ddt.f_p.eval           = function(X,Omega) apply(X,1,function(x) ddt.f_p(x,Omega))
ddx.ddt.f_p.eval       = function(X,Omega) apply(X,1,function(x) ddx.ddt.f_p(x,Omega))
ddOmega_p.ddt.f_p.eval = function(X,Omega) apply(X,1,function(x) ddOmega_p.ddt.f_p(x,Omega))


## logPost_p
# goal: log posterior density of the process model (first level of inference)
# logPost_p = function(Ybar_o,ddt.Ybar_o,Omega_p)
# {
#   ## predict
#   ddt.Ybar_p = ddt.f_p.eval(Ybar_o,Omega_p)
# 
#   ## residuals
#   res_p = ddt.Ybar_o - ddt.Ybar_p
# 
#   ## Likelihood
#   logLik_p   = - sum(res_p^2)/(sd1_p^2)
# 
#   ## prior
#   logPrior_p = - sum(Omega_p^2)/(sd2_p^2)
# 
#   ## posterior
#   logPost_p  = logLik_p + logPrior_p
# 
#   ## terminate
#   return(logPost_p)
# }
# 
# ## ddOmega_p.logPost_p
# # goal: derivative of the log posterior density of the process model wtr to parameters
# ddOmega_p.logPost_p = function(Ybar_o,ddt.Ybar_o,Omega_p)
# {
#   ## predict
#   ddt.Ybar_p = ddt.f_p.eval(Ybar_o,Omega_p)
# 
#   ## derivative of residuals
#   res_p           = ddt.Ybar_o - ddt.Ybar_p
#   ddOmega_p.res_p = - ddOmega_p.ddt.f_p.eval(Ybar_o,Omega_p)
# 
#   ## derivative of likelihood
#   ddOmega_p.logLik_p = - 2 * ddOmega_p.res_p%*%res_p/(sd1_p^2)
# 
#   ## derivative of prior
#   ddOmega_p.logPrior_p = - 2 * Omega_p /(sd1_p^2)
# 
#   ## derivative of posterior
#   ddOmega_p.logPost_p = ddOmega_p.logLik_p + ddOmega_p.logPrior_p
# 
#   ## terminate
#   return(ddOmega_p.logPost_p)
# }

## logPost_p
# goal: log posterior density of the process model (first level of inference)
logPost_p = function(Ybar_o,ddt.Ybar_o,Omega_p)
{
  ## predict
  ddt.Ybar_p = ddt.f_p.eval(Ybar_o,Omega_p)

  ## residuals
  res_p = ddt.Ybar_o - ddt.Ybar_p

  ## Likelihood
  logLik_p   = - 0.5 * length(ddt.Ybar_o) * log(0.5 * sum(res_p^2) + 1)

  ## prior
  logPrior_p = - 0.5 * length(Omega_p) * log(0.5 * sum(Omega_p^2) + 1)

  ## posterior
  logPost_p  = logLik_p + logPrior_p

  ## terminate
  return(logPost_p)
}

## ddOmega_p.logPost_p 
# goal: derivative of the log posterior density wtr to parameters
ddOmega_p.logPost_p = function(Ybar_o,ddt.Ybar_o,Omega_p)
{
	## predict
	ddt.Ybar_p = ddt.f_p.eval(Ybar_o,Omega_p)

	## derivative of residuals
	res_p = ddt.Ybar_o - ddt.Ybar_p
	ddOmega_p.res_p = -ddOmega_p.ddt.f_p.eval(Ybar_o,Omega_p)

	## derivative of likelihood
	ddOmega_p.logLik_p = - 0.5 * length(res_p) * 1/(0.5 * sum(res_p^2) + 1) * 0.5 * ddOmega_p.res_p%*%res_p

	## derivative of prior
	ddOmega_p.logPrior_p = - 0.5 * length(Omega_p) * 1/(0.5 * sum(Omega_p^2) + 1) * Omega_p

	## derivative of posterior
	ddOmega_p.logPost_p = ddOmega_p.logLik_p + ddOmega_p.logPrior_p

	## terminate
	return(ddOmega_p.logPost_p)
}

## argmax.logPost_p
# goal: find the maximum of the log posterior density of the process model
argmax.logPost_p = function(Ybar_o,ddt.Ybar_o,Omega_p)
{
	error_     = function(x) -logPost_p(Ybar_o,ddt.Ybar_o,x)
	graderror_ = function(x) -ddOmega_p.logPost_p(Ybar_o,ddt.Ybar_o,x)
	Omega_p  = optim(par=Omega_p,
					 fn=error_,
					 gr=graderror_,
					 method="BFGS"#,
					 # control=list("trace"=1,"REPORT"=1,"maxit"=100)
					 )$par
	return(Omega_p)
}

## .plot.DIN
## goal: plot the dynamical interaction network of the system
# effectsMat - matrix - matrix of pairwise effects between system variables (e.g. row 1 col 2 is the effect of variable 2 on variable 1)
# weightsMat - matrix - matrix of pairwise weights of the effects between system variables (e.g. row 1 col 2 corresponds to the contribution of variable 2 on variable 1)
# labels     - vector - vector of the names of the variables in the matrix
.plot.DIN = function(effectsMat,weightsMat,labels)
{
	## dimensions
	N = dim(effectsMat)[1]

	## scale effects and contributions
	effectsMat = (effectsMat>0)*1
	# weightsMat = weightsMat/sum(weightsMat) # proportion of total change

	## angles
	theta = seq(0,2*pi,(2*pi)/N)
	x = cos(theta)
	y = sin(theta)
	x_ = cos(theta+(2*pi*0.05))
	y_ = sin(theta+(2*pi*0.05))
	x__ = 1.25*cos(theta+(2*pi*0.025))
	y__ = 1.25*sin(theta+(2*pi*0.025))

	## plot interactions
	plot(x=c(-1:1)*2,y=c(-1:1)*2,cex=0,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
	for(i in 1:N)
	{
		for(j in 1:N)
		{
			color_ = if(effectsMat[i,j]>0){"green"}else{"red"}
			# points(x__[i],y__[i],cex=30/N)
			text(x__[i],y__[i],labels=labels[i])
			if(weightsMat[i,j]*10>0)
			{
				arrows(x0=x[j],x1=x_[i],y0=y[j],y1=y_[i],lwd=weightsMat[i,j]*10,col=color_,length=0.1)
			}
		}
	}
}

#
###

###############
## FIT MODEL ##
###############

## goal: fit the model to the time series

## initiate model
n     = nrow(TS)
n_    = n*4
N     = ncol(TS) - 1
N_e   = 30
W_o   = 10
W_p   = 10
sd1_o = 0.01 # only relevant for first level of inference
sd2_o = .4 
sd1_p = 0.01 # only relevant for first level of inference
sd2_p = 1.0 
alpha = 10
t     = TS[,1]
Y     = TS[,-1]
t     = (t-min(t))/(max(t)-min(t))*alpha
mean_ = apply(Y,2,mean)
sd_   = apply(Y,2,sd)
Y     = t(apply(Y,1,function(x)(x-mean_)/sd_)) # standardise Y
t_    = seq(min(t),max(t),(max(t)-min(t))/n_)
n_    = length(t_)

## fit observation model
i = 1
Ybar_o     = NULL
ddt.Ybar_o = NULL
for(i in 1:N)
{
	## get ensemble
    ensemble_ = NULL
    for(e in 1:N_e)
    {
      Omega_o     = rnorm(3*W_o,0,sd2_o)
      Omega_o     = argmax.logPost_o(t,Y[,i],Omega_o)
      ensemble_   = rbind(ensemble_,Omega_o)
    }

    ## get Ybar
    f_o.eval_ = apply(ensemble_,1,function(x)f_o.eval(t_,x))
    Ybar_o_q0.05 = apply(f_o.eval_,1,quantile,probs=0.05)
    Ybar_o_q0.25 = apply(f_o.eval_,1,quantile,probs=0.25)
    Ybar_o_q0.50 = apply(f_o.eval_,1,quantile,probs=0.50)
    Ybar_o_q0.75 = apply(f_o.eval_,1,quantile,probs=0.75)
    Ybar_o_q0.95 = apply(f_o.eval_,1,quantile,probs=0.95)
    Ybar_o_mean  = apply(f_o.eval_,1,mean)

    ## get ddt.Ybar
    ddt.f_o.eval_ = apply(ensemble_,1,function(x)ddt.f_o.eval(t_,x))
    ddt.Ybar_o_q0.05 = apply(ddt.f_o.eval_,1,quantile,probs=0.05)
    ddt.Ybar_o_q0.25 = apply(ddt.f_o.eval_,1,quantile,probs=0.25)
    ddt.Ybar_o_q0.50 = apply(ddt.f_o.eval_,1,quantile,probs=0.50)
    ddt.Ybar_o_q0.75 = apply(ddt.f_o.eval_,1,quantile,probs=0.75)
    ddt.Ybar_o_q0.95 = apply(ddt.f_o.eval_,1,quantile,probs=0.95)
    ddt.Ybar_o_mean  = apply(ddt.f_o.eval_,1,mean)

    ## visualise fit
    par(mfrow=c(2,1))
    #
    plot(t,Y[,i],type="b",ylim=c(-3,3))
    polygon(c(t_,rev(t_)),c(Ybar_o_q0.05,rev(Ybar_o_q0.95)),col=rainbow(N,alpha=0.2)[i],border=NA)
    polygon(c(t_,rev(t_)),c(Ybar_o_q0.25,rev(Ybar_o_q0.75)),col=rainbow(N,alpha=0.2)[i],border=NA)
    lines(t_,Ybar_o_q0.50,col=rainbow(N)[i],lty=2)
    lines(t_,Ybar_o_mean,col=rainbow(N)[i],lty=3)
    for(k in 1:ncol(f_o.eval_))
    {
      lines(t_,f_o.eval_[,k],col=rainbow(N,alpha=0.2)[i])
    }
    #
    plot(t_,ddt.Ybar_o_q0.50,ylim=c(-3,3),cex=0)
    polygon(c(t_,rev(t_)),c(ddt.Ybar_o_q0.05,rev(ddt.Ybar_o_q0.95)),col=rainbow(N,alpha=0.2)[i],border=NA)
    polygon(c(t_,rev(t_)),c(ddt.Ybar_o_q0.25,rev(ddt.Ybar_o_q0.75)),col=rainbow(N,alpha=0.2)[i],border=NA)
    lines(t_,ddt.Ybar_o_q0.50,col=rainbow(N)[i],lty=2)
    lines(t_,ddt.Ybar_o_mean,col=rainbow(N)[i],lty=3)
    for(k in 1:ncol(ddt.f_o.eval_))
    {
      lines(t_,ddt.f_o.eval_[,k],col=rainbow(N,alpha=0.2)[i])
    }
    #
    par(mfrow=c(1,1))

    ## store
    Ybar_o = cbind(Ybar_o,Ybar_o_mean)
    ddt.Ybar_o = cbind(ddt.Ybar_o,ddt.Ybar_o_mean)
}

## fit process model
EMat = NULL
CMat = NULL
r2Vect   = NULL
ensemble = NULL
for(i in 1:N)
{
  ## BFGS
  ensemble_ = NULL
  for(e in 1:N_e)
  {
    Omega_p   = rnorm((2+N)*W_p,0,sd2_p)
    Omega_p   = argmax.logPost_p(Ybar_o,ddt.Ybar_o[,i],Omega_p)
    ensemble_ = rbind(ensemble_,Omega_p)
  }
  ensemble = cbind(ensemble,ensemble_)
}

## analyse process model
# pdf(paste("TS_analysis_",m,".pdf",sep=""),height=18,width=6)
#
for(i in 1:N)
{
  ## get ensemble for a specific time series
  D = (2+N)*W_p
  ensemble_ = ensemble[,D*(i-1)+1:D]
  
  ## compute ddt.Ybar_p
  tmp = apply(ensemble_,1, function(x) ddt.f_p.eval(Ybar_o,x))
  ddt.Ybar_p_q0.05 = apply(tmp,1,quantile,probs=0.05)
  ddt.Ybar_p_q0.50 = apply(tmp,1,quantile,probs=0.50)
  ddt.Ybar_p_q0.95 = apply(tmp,1,quantile,probs=0.95)
  ddt.Ybar_p_mean = apply(tmp,1,mean)

  ## compute effects
  tmp = apply(ensemble_,1, function(x) ddx.ddt.f_p.eval(Ybar_o,x))
  ddx.ddt.Ybar_p_q0.05 = matrix(apply(tmp,1,quantile,probs=0.05),nrow=n_,byrow=T)
  ddx.ddt.Ybar_p_q0.50 = matrix(apply(tmp,1,quantile,probs=0.50),nrow=n_,byrow=T)
  ddx.ddt.Ybar_p_q0.95 = matrix(apply(tmp,1,quantile,probs=0.95),nrow=n_,byrow=T)
  ddx.ddt.Ybar_p_mean = matrix(apply(tmp,1,mean),nrow=n_,byrow=T)

  ## compute mean effects
  Emean_q0.50 = apply(ddx.ddt.Ybar_p_q0.50,2,mean)
  Emean_mean = apply(ddx.ddt.Ybar_p_mean,2,mean)

  ## compute contributions
  C_q0.05 = ddx.ddt.Ybar_p_q0.05*ddt.Ybar_o
  C_q0.50 = ddx.ddt.Ybar_p_q0.50*ddt.Ybar_o
  C_q0.95 = ddx.ddt.Ybar_p_q0.95*ddt.Ybar_o
  C_mean = ddx.ddt.Ybar_p_mean*ddt.Ybar_o

  ## compute total contributions
  Ctot_q0.05 = t(apply(C_q0.05^2,1,function(x)x/sum(x)))
  Ctot_q0.50 = t(apply(C_q0.50^2,1,function(x)x/sum(x)))
  Ctot_q0.95 = t(apply(C_q0.95^2,1,function(x)x/sum(x)))
  Ctot_mean = t(apply(C_mean^2,1,function(x)x/sum(x)))
  

  ## r squared
  res = ddt.Ybar_o[,i] - ddt.Ybar_p_mean
  r2  = 1 - (sd(res)^2)/(sd(ddt.Ybar_o[,i])^2)

  EMat = rbind(EMat,Emean_mean)
  CMat = rbind(CMat,apply(Ctot_mean,2,sum))
  r2Vect = c(r2Vect,r2)


  par(mfrow=c(5,1),mar=c(4,5,1,1),cex.lab=2)
  #
  ## plot time series
  plot(t_,Ybar_o[,i],type="l",col=rainbow(N)[i],xlab="Time",ylab="Density")
  points(t,Y[,i])
  #
  ## plot dynamics
  plot(t_,ddt.Ybar_o[,i],type="l",xlab="Time",ylab="Dynamics")
  polygon(c(t_,rev(t_)),c(ddt.Ybar_p_q0.05,rev(ddt.Ybar_p_q0.95)),border=NA,col=rainbow(N,alpha=0.2)[i])
  lines(t_,ddt.Ybar_p_mean,type="l",col=rainbow(N)[i])
  lines(t_,rep(0,n_),lty=2)
  #
  ## plot effects
  plot(t_,ddx.ddt.Ybar_p_mean[,1],ylim=c(-5,5),cex=0,xlab="Time",ylab="Effect")
  for(j in 1:N)
  {
    polygon(c(t_,rev(t_)),c(ddx.ddt.Ybar_p_q0.05[,j],rev(ddx.ddt.Ybar_p_q0.95[,j])),border=NA,col=rainbow(N,alpha=0.2)[j])
    lines(t_,ddx.ddt.Ybar_p_mean[,j],type="l",col=rainbow(N)[j])
  }
  lines(t_,rep(0,n_),lty=2)
  #
  ## plot contributions
  plot(t_,C_mean[,1],ylim=c(-6,6),cex=0,xlab="Time",ylab="Contribution")
  for(j in 1:N)
  {
    polygon(c(t_,rev(t_)),c(C_q0.05[,j],rev(C_q0.95[,j])),border=NA,col=rainbow(N,alpha=0.2)[j])
    lines(t_,C_mean[,j],type="l",col=rainbow(N)[j])
  }
  lines(t_,rep(0,n_),lty=2)
  #
  ## plot contributions (%)
  X_q0.05 = Ctot_q0.05
  X_q0.50 = Ctot_q0.50
  X_q0.95 = Ctot_q0.95
  X_mean = Ctot_mean
  plot(t_,X_q0.50[,1],ylim=c(0,1),cex=0,xlab="Time",ylab="Contribution (%)")
  for(j in 1:N)
  {
    # polygon(c(t_,rev(t_)),c(X_q0.05[,j],rev(X_q0.95[,j])),border=NA,col=rainbow(N,alpha=0.2)[j])
    lines(t_,X_mean[,j],type="l",col=rainbow(N)[j])
  }
  lines(t_,rep(0,n_),lty=2)
  #
  par(mfrow=c(1,1),mar=c(5,4,3,3)+0.1)
  
}
#
# dev.off()

# ## compute equilibrium densities
# argmin.ddt.f_p = function(Y_0,Omega)
# {
#   ## minimise ddt.f_p
#   Omega = matrix(Omega,ncol=N,byrow=F)
#   fn = function(x) ddt.f_p(x, Omega[,1])^2 + ddt.f_p(x, Omega[,2])^2 + ddt.f_p(x, Omega[,3])^2
#   res = optim(par=Y_0, fn=fn)
#   return(c(res$value,res$par))
# }
# #
# tmp = apply(ensemble,1, function(x) argmin.ddt.f_p(rep(0,N),x))
# tmp = tmp[,which(tmp[1,]<=0.001)][-1,]
# for(i in 1:N) 
# {
#    	torm = which(abs(tmp[i,])>=3) 
# 	if( length(torm) != 0 )
# 	{
# 		tmp = tmp[,-torm]
# 	}
# }
# argmin.ddt.f_p_q0.05 = matrix(apply(tmp,1,quantile,probs=0.05),nrow=N,byrow=T)
# argmin.ddt.f_p_q0.50 = matrix(apply(tmp,1,quantile,probs=0.50),nrow=N,byrow=T)
# argmin.ddt.f_p_q0.95 = matrix(apply(tmp,1,quantile,probs=0.95),nrow=N,byrow=T)
# #
# tmp = t(tmp); colnames(tmp) = colnames(Y)
# boxplot(tmp,ylab="Density")

##
# pdf(paste("TS_DIN_",m,".pdf",sep=""),height=6,width=6)
#
CMat = t(apply(CMat,1,function(x)x/sum(x)))
.plot.DIN(EMat,CMat*(CMat>0.1),colnames(Y))
#
# dev.off()

## terminate
outputList = list(r2Vect,EMat,CMat)

#
###
