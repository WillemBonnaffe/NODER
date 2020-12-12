#############
## NODER.R ##
#############

## goal: declare all functions involved in the NODE fitting process

## update
## 12-12-2020 - created NODER_v1.0.0.r


##############
## INITIATE ##
##############

## load ODE library
library(deSolve)

#
###

################
## TIME SERIE ##
################

## goal: functions to handle time series data

## TS.format
## goal: format a time series by transforming then standardising variables 
# TS       - matrix - matrix containing system variables including time by columns
# indexLog - vector - vector containing index of variables to log 
TS.format = function(TS,indexLog=NULL)
{
	## standardise years
	TS[,1] = TS[,1] - min(TS[,1])
	# TS[,1] = (TS[,1] - mean(TS[,1]))/sd(TS[,1]) # to standardise year
	
	## transform positive variables
	if(!is.null(indexLog))
	{
		TS[,indexLog] = log(TS[,indexLog])
	}
		
	## standardise
	mean_ = apply(TS[,-1],2,mean,na.rm=T)                       
	sd_ = apply(TS[,-1],2,sd,na.rm=T)                           
	TS[,-1] = t(apply(TS[,-1],1,function(x){(x - mean_)/sd_}))

	## terminate
	return(TS)
}

## TS.differentiate
## goal: compute time difference between consecutive measurements in each variable of the time series
# TS - matrix - matrix containing system variables including time by columns
TS.timedifferentiate = function(TS)
{
	## initiate
	t = TS[,1]
	Y = TS[,-1]
	N = nrow(TS)

	## time difference 
	dt = unique(round(diff(t),4))
	if(length(dt)!=1){message("warning: dt not constant")}

	## compute time difference
	dY = apply(Y,2,diff)
	dY = dY/dt 

	## combine
	dTS = cbind(t[-N],dY)
	colnames(dTS) = colnames(TS)

	## terminate
	return(dTS)
}

#
###

#############################
## SINGLE LAYER PERCEPTRON ##
#############################

## goal: functions to implement single layer perceptrons in R (i.e. 1 hidden layer neural networks)

## list of activation function types for SLPs
f_sigmaList = list(
	f_sigma = function(X){return(X)},
	f_sigma = function(X){return(1/(1+exp(-X)))},
	f_sigma = function(X){return(exp(-X^2))},
	f_sigma = function(X){return((X>0)*X)},
	f_sigma = function(X){return(exp(X))},
	f_sigma = function(X){return(sin(X))}
)

## list of derivatives of activation functions for SLPs
df_sigmaList = list(
	df_sigma = function(X){return(1)},
	df_sigma = function(x){return((1-(1/(1+exp(-x))))*(1/(1+exp(-x))))},
	df_sigma = function(X){return(-2*X*exp(-X^2))},
	df_sigma = function(X){return(X>0)},
	df_sigma = function(X){return(exp(X))},
	df_sigma = function(X){return(cos(X))}
)

## SLP 
## goal: evaluate a SLP for given arguments
# X       - vector - vector of inputs/covariates
# Omega   - vector - vector of parameters of the SLP
# spec    - vector - vector containing [1] input size, [2] width of the hidden layer (i.e. number of hidden nodes), [3] output size
# f_sigma - func   - activation function
SLP = function(X,Omega,spec,f_sigma)
{
	## dimensions
	I = spec[1]
	W = spec[2]
	O = spec[3]
	
	## weight vectors
	w0 = Omega[1:O]
	w1 = t(matrix(Omega[1:(O*W)+O],ncol=O))
	w2 = Omega[1:W+O*W+O]
	w3 = matrix(Omega[1:(I*W)+W+O*W+O],ncol=I)

	## evaluate SLP
	res = w0 + w1%*%f_sigma(w2 + w3%*%X)
	return(res)
}

## dSLP
## goal: evaluate the derivative of the SLP
# X        - vector - vector of inputs/covariates
# Omega    - vector - vector of parameters of the SLP
# spec     - vector - vector containing [1] input size, [2] width of the hidden layer (i.e. number of hidden nodes), [3] output size
# df_sigma - func   - derivative of the activation function
dSLP = function(X,Omega,spec,df_sigma)
{	
	## dimensions
	I = spec[1]
	W = spec[2]
	O = spec[3]

	## weight vectors
	w0 = Omega[1:O]
	w1 = t(matrix(Omega[1:(O*W)+O],ncol=O))
	w2 = Omega[1:W+O*W+O]
	w3 = matrix(Omega[1:(I*W)+W+O*W+O],ncol=I)

	## evaluate gradient	
	dSLP_ = w1%*%(w3*as.vector(df_sigma(w2 + w3%*%X)))
	return(dSLP_)
}

## SLP.dim
## goal: get the total number of parameters in the SLP
# spec - vector - vector containing [1] input size, [2] width of the hidden layer (i.e. number of hidden nodes), [3] output size
SLP.dim = function(spec)
{
	## dimensions
	I = spec[1]
	W = spec[2]
	O = spec[3]
	return(I*W+W+O*W+O)
}

#
###

####################
## NODE FUNCTIONS ##
####################

## goal: functions to define and manipulate a NODE model

## NODE.init
## goal: create a NODE object list based on user-defined NODE specifications 
# TS                - matrix - matrix containing the time series to be analysed with the NODE model, the first column should be time
# SLPinputStateList - list   - list containing the index of the system state variables that each equation in the NODE system depends on
# SLPftypeVect      - vector - vector containing the activation function type for each NODE
NODE.init = function(TS,SLPinputStateList,SLPftypeVect)
{
	## variables
	# t         - vector - vector of time steps over which to evaluate the associated function
	# Y_0       - vector - vector of initial system state variable values excluding time
	# state     - vector - vector of system variable states including time
	# Omega     - vector - is the parameter vector for all NODE
	# OmegaList - list   - list contaiing the parameter vectors for each NODE 
	# Theta     - vector - is the overall parameter vector containing initial conditions and Omega, i.e. Theta = [Y_0,Omega]

	################################
	## -- SYSTEM SPECIFICATION -- ##

	## goal: specify NODE system (dYdt equations) and Jacobian

	## system properties
	t      = TS[,1]       # times
	Y      = TS[,-1]      # system state variables excluding time 
	n      = nrow(TS)     # number of rows
	N      = ncol(TS)     # number of system state variables including time
	i_init = 1:(N-1)      # index of initial values in main parameter vector Theta

	## initiate SLP specifications (n input, n nodes, n outputs) and dimension (n parameters)
	# i = 1 is time => SLP attributes of index 1 are irrelevant => set to NA
	SLPinputStateList = c(NA,SLPinputStateList)  # list containing the index of the system state variables that each equation in the NODE system depends on
	SLPftypeVect      = c(NA,SLPftypeVect)       # vector specifying the type of activation function to be used for each NODE (e.g. 1 is identity, 2 is sigmoid, ..., see SLP function section)
	SLPspecMat        = matrix(NA,nrow=N,ncol=3) # matrix containing the specficiation of each SLP in the NODEs listed by row, each row contains (1) number inputs, (2) number hidden units, (3) number of outputs
	SLPdimVect        = 0                        # vector containing the number of parameters in each SLP
	SLPfsigmaList     = list(NA)                 # list containig the activation functions
	SLPdfsigmaList    = list(NA)                 # list containgin the derivative of the activation functions
	W = 10                                       # number of hidden units in all SLPs
	O = 1                                        # number of outputs in all SLPs

	## compute specifications for each NODE
	for(i in 2:N) 	
	{
		I                   = length(SLPinputStateList[[i]]) # number of inputs
		SLPspecMat[i,]      = c(I,W,O)
		SLPdimVect[i]       = SLP.dim(c(I,W,O))
		SLPfsigmaList[[i]]  = f_sigmaList[[SLPftypeVect[i]]]
		SLPdfsigmaList[[i]] = df_sigmaList[[SLPftypeVect[i]]]
	}

	## NODE.Ydot
	## goal: NODE equation system with each NODE being defined as a SLP, i.e. Ydot_i = SLP_i(state,Omega)
	NODE.Ydot = function(state,OmegaList)
	{	
		## initiate 
		Ydot = 1 # i = 1 is time => dY/dt[1] = 1

		## for each NODE
		for(i in 2:N) 
		{
			Ydot[i] = SLP(state[SLPinputStateList[[i]]],OmegaList[[i]],SLPspecMat[i,],SLPfsigmaList[[i]])
		}
		return(Ydot)
	}

	## NODE.J
	## goal: compute Jacobian matrix associated withe the NODE system above, i.e. J = [d(dY_i/dt)/dY_j]_ij 
	NODE.J = function(state,OmegaList)
	{	
		## initiate
		J = matrix(0,N,N) # i = 1 is time => J[1,] = d(dY/dt_1)/dY_1, d(dY/dt_1)/dY[2], ... = 0

		## for each NODE
		for(i in 2:N) 
		{
			stateIndexVect      = SLPinputStateList[[i]]
			J[i,stateIndexVect] = dSLP(state[stateIndexVect],OmegaList[[i]],SLPspecMat[i,],SLPdfsigmaList[[i]])
		}
		return(J)
	}

	## compute number of parameters for each NODE
	dimVect = SLPdimVect         # the number of parameters in each DE depends on the DE structure specified
	D       = sum(dimVect) + N-1 # overall number of parameters in the NODE equations + initial conditions
	
	################################
	## -- ENGINE SPECIFICATION -- ##

    ## goal: specify simulation engine

	## NODE.Ybar.tde
	## goal: tde (time difference equation) simulation engine
	NODE.Ybar.tde = function(t,Y_0,OmegaList)
	{
		Ybar = .eval.tde(t,Y_0,function(x)NODE.Ydot(x,OmegaList)[-1])
		return(Ybar)
	}

	## NODE.Ybar.ode
	## goal: ode simulation engine (ordinary differential equation)
	NODE.Ybar.ode = function(t,Y_0,OmegaList)
	{
		Ybar = .eval.ode(t,Y_0,function(x)NODE.Ydot(x,OmegaList)[-1])
		return(Ybar)
	}

	######################################
	## FUNCTIONS FOR MANIPULATING MODEL ##

	## goal: define user friendly wrapers for the NODE. functions

	## general functions dYdt and Jacobian
	.s           = function(Omega)        .slice(Omega,dimVect)
	.Ydot        = function(state,Omega)  NODE.Ydot       (state, .s(Omega))
	.J           = function(state,Omega)  NODE.J          (state, .s(Omega))
	.C           = function(state,Omega)  t(t(NODE.J(state, .s(Omega)))*NODE.Ydot(state, .s(Omega)))

	## essential functions for time difference equations
	.Ybar.tde      = function(t,Y_0,Omega)  NODE.Ybar.tde    (t,Y_0, .s(Omega))
	.Ydotbar.tde   = function(t,Y_0,Omega)  t(apply(.Ybar.tde(t,Y_0,Omega),1,function(state).Ydot(state,Omega)))
	.Jbar.tde      = function(t,Y_0,Omega)  t(apply(.Ybar.tde(t,Y_0,Omega),1,function(state).J(state,Omega)))
	.Cbar.tde      = function(t,Y_0,Omega)  t(apply(.Ybar.tde(t,Y_0,Omega),1,function(state).C(state,Omega)))

	## essential functions for ordinary differential equations
	.Ybar.ode      = function(t,Y_0,Omega)  NODE.Ybar.ode   (t,Y_0, .s(Omega))
	.Ydotbar.ode   = function(t,Y_0,Omega)  t(apply(.Ybar.ode(t,Y_0,Omega),1,function(state).Ydot(state,Omega)))
	.Jbar.ode      = function(t,Y_0,Omega)  t(apply(.Ybar.ode(t,Y_0,Omega),1,function(state).J(state,Omega)))
	.Cbar.ode      = function(t,Y_0,Omega)  t(apply(.Ybar.ode(t,Y_0,Omega),1,function(state).C(state,Omega)))

	###############################
	## -- INSTANTIATE OBJECTS -- ##

	## goal: create NODE object list

	## initiate parameter
	Theta = rnorm(D,0,0.01)

	## create NODE model 
	model = list()

	## properties
	model$TS                    = TS                   # matrix - the time series with time in first column
	model$n                     = n                    # int    - number of rows in the time series
	model$N                     = N                    # int    - number of columns in the time series (i.e. number of variables in the NODE system)
	model$i_init                = i_init               # vector - indices of initial state values (excluding time) in the main parameter vector, it should be the N-1 first elements
	model$SLPinputStateList     = SLPinputStateList    # list   - each element is a vector specifying the dependence of dYdt_i on other variables j in the system
	model$SLPftypeVect          = SLPftypeVect         # vector - type of activation function to use for each NODE in the system
	model$SLPfsigmaList         = SLPfsigmaList        # list   - activation function of each NODE in the system (the first element is time)
	model$SLPdfsigmaList        = SLPdfsigmaList       # list   - derivative of activation function of each NODE in the system
	model$SLPspecMat            = SLPspecMat           # matrix - 3 columns, one for number of inputs, number of hidden nodes, number of outputs
	model$dimVect               = dimVect              # vector - number of parameters in each NODE
	model$D                     = D                    # int    - total number of parameters in the NODE model (including initial state values)
	model$Theta                 = Theta                # vector - main parameter vector containing all parameters in the NODE

	## general 
	model$.s                    = .s                   # func   - slice function to distribute the parameter vector into a list of parameters for each of the NODE
	model$.Ydot                 = .Ydot                # func   - returns temporal change in each state variable (i.e. dY_i/dt), takes state vector and equation parameter vector (Omega) as input
	model$.J                    = .J                   # func   - returns the Jacobian of the NODE system
	model$.C                    = .C                   # func   - returns the product of the Jacobian and the dynamics of each state variebles (i.e. Geber method)

	## essential tde
	model$.Ybar.tde             = .Ybar.tde            # func   - simulate NODE system and returns time series state variables (uses a time difference simulation scheme)
	model$.Ydotbar.tde          = .Ydotbar.tde         # func   - evaluate NODE system (i.e. dYdt) at every time step (uses a time difference simulation scheme)
	model$.Jbar.tde             = .Jbar.tde            # func   - evaluate Jacobian at every time step (uses a time difference simulation scheme)
	model$.Cbar.tde             = .Cbar.tde            # func   - evaluate contributions at every time step (uses a time difference simulation scheme)

	## essential ode
	model$.Ybar.ode             = .Ybar.ode            # func   - simulate NODE system and returns time series state variables (uses an ode simulation scheme)
	model$.Ydotbar.ode          = .Ydotbar.ode         # func   - evaluate NODE system (i.e. dYdt) at every time step (uses an ode simulation scheme)
	model$.Jbar.ode             = .Jbar.ode            # func   - evaluate Jacobian at every time step (uses an ode simulation scheme)
	model$.Cbar.ode             = .Cbar.ode            # func   - evaluate contributions at every time step (uses an ode simulation scheme)

	## terminate
	return(model)
}

## .eval.tde
## goal: simple time difference engine for simulating an ODE system 
# t    - vector - time steps to simulate the ODE over
# Y_0  - vector - initial value of system state variables 
# dYdt - func   - temporal change in each state variable (i.e. dY_i/dt), input should be the state variable vector (including time as first element) 
.eval.tde = function(t,Y_0,dYdt)
{
	## initiate 
	n        = length(t)
	N        = length(Y_0)
	Ybar     = matrix(NA,nrow=n,ncol=N)
	Ybar[1,] = Y_0

	## euler step scheme
	for(k in 2:n)
	{
		t_       = t[(k-1):k]
		Ybar_    = Ybar[k-1,]
		Ybar[k,] = Ybar_ + (t_[2]-t_[1])*dYdt(c(t_[1],Ybar_))
	}

	## terminate
	return(cbind(t,Ybar))
}

## .eval.ode
## goal: numerical ode solver engine for simulating an ODE system 
# t    - vector - time steps to simulate the ODE over
# Y_0  - vector - initial value of system state variables 
# dYdt - func   - temporal change in each state variable (i.e. dY_i/dt), input should be the state variable vector (including time as first element) 
.eval.ode = function(t,Y_0,dYdt)
{
	## odeSolve
	fnwrapper = function(time,y,parms){list(dYdt(c(time,y)))}
	return(ode(y=Y_0, times=t, func=fnwrapper, parms=NULL, method="ode45", maxsteps=100))
}

## .error.LR
## goal: 
# Y      - matrix - response matrix with predicted variables by columns
# X      - matrix - predictor matrix with explanatory variables by columns
# f      - func   - predictive function to map X on Y, i.e. Ybar = f(X), where Ybar is the predicted value of the response
# param  - vector - vector of parameter values controlling the predictive function
# lambda - float  - scaling factor controlling the penalisation of the coefficient values (high values correspond to high penalisation)
# alpha  - prop   - parameter between 0 and 1 to control the mixture between LASSO and ridge penalization (0 for LASSO, 1 for ridge, 0.5 for 50% mixture)
.error.LR = function(Y,X,f,param,lambda,alpha=1)
{
	res     = as.vector(t(Y - f(X,param)))
	MSE     = mean((res)^2)
	penalty = alpha*lambda*sum(abs(param)) + (1-alpha)*lambda*sum(param^2)
	error   = MSE + penalty
	return(list(error,MSE,penalty))
}

## .error.normalBayesian
## goal: 
# Y        - matrix - response matrix with predicted variables by columns
# X        - matrix - predictor matrix with explanatory variables by columns
# f        - func   - predictive function to map X on Y, i.e. Ybar = f(X), where Ybar is the predicted value of the response
# param    - vector - vector of parameter values controlling the predictive function
# sd_lik   - vector - vector of standard deviation in the likelihood function for each response variable in the Y matrix
# mu_prior - vector - means of the normal prior distribution on each parameter
# mu_prior - vector - standard deviation of the normal prior distribution on each parameter
.error.normalBayesian = function(Y,X,f,param,sd_lik,mu_prior,sd_prior)
{
	res = as.vector(t(Y - f(X,param)))
	logLik   = sum(log(dnorm(res,0,sd_lik)))
	logPrior = sum(log(dnorm(param,mu_prior,sd_prior)))
	logPost  = logLik + logPrior
	return(list(logPost,logLik,logPrior))
}

## .fit
## goal: fit NODE model with BFGS gradient descent
# paramVect - vector - parameter values to fit
# error     - func   - function to evaluate error associated with a given parameter vector, takes the parameter vector as single input
# nIt       - int    - number of iterations of the optimisation
.fit = function(paramVect,error,nIt)
{
	## run BFGS
	optimised = optim(par=paramVect,fn=error,method="BFGS",control=list("trace"=1, "REPORT"=1, "maxit"=nIt)) 
	return(optimised)
}

## .ABC
## goal: approximate posterior distribution with rejection sampling
# paramVect_0 - vector - vector of initial parameter values to perturbate
# error       - func   - function to evaluate error essociated with a given parameter vector, takes the parameter vector as single input
# nIt         - int    - number of iterations of the optimisation
# noise       - float  - standard distance of the uniform perturbation applied on the initial parameter vector 
# threshold   - prop   - proportion of the samples that should be kept (e.g. 0.5 means 50% sample rejected) 
.ABC = function(paramVect_0,error,nIt,noise,threshold)
{
	## sample
	chain = NULL
	for(k in 1:nIt)
	{
		## uniform perturbation on parameters
		paramVect_ = paramVect_0 + runif(length(paramVect_0),0,noise)

		## compute error
		error_ = error(paramVect_)

		## store in chain
		chain = rbind(chain,c(error_,paramVect_))
	}

	## reject samples x%
	threshold = quantile(chain[,1],probs=threshold)
	chain     = chain[which(chain[,1]>=threshold),]

	## terminate 
	chain = rbind(c(error(paramVect_0),paramVect_0),chain) 
	MaP   = chain[which.max(chain[,1]),-1]
	return(list("MaP"=MaP,"chain"=chain))
}

## .plot.fit
## goal: visualise a time series object along with a predicted time series object (e.g. from an ODE model) 
# TS         - matrix - time series matrix containing time in first column and system state variables in the other columns
# prediction - matrix - time series matrix containing time in first column and predicted system state variables in other columns
.plot.fit = function(TS,prediction=NULL)
{
	## intiate
	N = ncol(TS)

	## plot pairwise variables
	colVect = c("black",rainbow(N-1))
	par(mfrow=c(N,N),oma=c(2,2,1,1),mar=c(1,1,1,1)*0.5)
	for(i in 1:N)
	{
		for(j in 1:N)
		{
			if(i>j)
			{
				## set axes 
				if((j==1)&(i<N)){xaxt="n";yaxt="s"}     # first plot column
				if((i==N)&(j>1)){xaxt="s";yaxt="n"}     # last plot row
				if((i==N)&(j==1)){xaxt="s";yaxt="s"}    # bottom left corner
				
				## plot
				plot(TS[,j],TS[,i],type="b",pch=16,xaxt=xaxt,yaxt=yaxt,col=adjustcolor("black",alpha=0.5))

				## plot best fit
				if(!is.null(prediction)) lines(prediction[,j],prediction[,i],col="red")
				
			} else
			{
				if(i == j)
				{
					plot(-1:1,-1:1,cex=0,xlab="",ylab="",xaxt="n",yaxt="n")
					text(0,0,colnames(TS)[i],cex=2)
				} else
				{
					plot.new()
				}
			}
		}

	}	
	#
	par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(5,4,4,2)+0.1) # default
}

## .residualCheck1.plot
## goal: visualise residuals for pairwise combinations of variables
# res - vector - a vector of residuals (i.e. differences between observed and predicted states)
.plot.residualCheck1 = function(res)
{
	## correlations between residuals
	plot(res)
	return(cor(res))
}

## .residualCheck2.plot
## goal: visualise autocorrelation in residuals for each system variable
# res - vector - a vector of residuals (i.e. differences between observed and predicted states)
.plot.residualCheck2 = function(res)
{
	## autocorrelations between residuals
	par(mfrow=c(ncol(res),1))
	#
	for(i in 1:ncol(res))
	{
		barplot(res[,i],main=colnames(res)[i])
	}
	#
	par(mfrow=c(1,1))
}

## .residualCheck3.plot
## goal: visualise standard qqplot on model residuals of each variable 
# res - vector - a vector of residuals (i.e. differences between observed and predicted states)
.plot.residualCheck3 = function(res)
{
	## qqplot
	par(mfrow=c(ncol(res),1))
	#
	for(i in 1:ncol(res))
	{
		q_theoretical = qnorm(p=seq(0,1,0.01),mean=mean(res[,i]),sd=sd(res[,i]))
		q_observed = quantile(x=res[,i],probs=seq(0,1,0.01))
		plot(q_observed,q_theoretical,main=colnames(res)[i])
		lines(c(-2,2),c(-2,2))
	}
	#
	par(mfrow=c(1,1))
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
	weightsMat = weightsMat/sum(weightsMat) # proportion of total change

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
			points(x__[i],y__[i],cex=30/N)
			text(x__[i],y__[i],labels=labels[i])
			arrows(x0=x[j],x1=x_[i],y0=y[j],y1=y_[i],lwd=weightsMat[i,j]*10,col=color_,length=0.1)
		}
	}
}

## .plot.PPS
## goal: plot a pairwise phase space of specified system state vaiable of the model for a single pairwise combination of variables
# xIndex - int    - index of the state variable to be plotted on the x axis
# yIndex - int    - index of the state variable to be plotted on the y axis
# zIndex - int    - index of the state variables to be considered as the response
# dYdt   - func   - function for the dyanmics of the state variables, takes a state vector containing only the two x and y states to plot over
# lims   - matrix - matrix containing the upper and lower limits of the x and y variables listed by rows
.plot.PPS = function(xIndex,yIndex,zIndex,dYdt,lims,labels=c("",""))
{

	## initiate
	MatList = list()
	colVect = rainbow(100,start=0,end=0.6,alpha=0.25)
	
	## for each variable in the system but time
	for(z in zIndex)
	{
		dYdt_ = function(x){return(dYdt(x)[z])}
		MatList[[z-1]] = .getMat(lims,dYdt_)
	}

	## discretise matrices into 0/1 (negative,positive)
	MatListGrey = MatList
	MatList = lapply(MatList,function(x)(x>0)*1)

	## overlay the matrices
	Mat = MatList[[1]]
	MatGrey = MatListGrey[[1]]
	for(i in 1:length(MatList))
	{
		MatGrey = MatGrey*MatListGrey[[i]]
		Mat = Mat + MatList[[i]]
	}

	## heatmap
	# image(MatGrey,col=adjustcolor(colVect,alpha=0.5))
	image(Mat,col=colVect,xaxt="n",yaxt="n",xlab=labels[1],ylab=labels[2])
	axis(1,labels=seq(min(lims[1,]),max(lims[1,]),(max(lims[1,])-min(lims[1,]))/4),at=seq(0,1,0.25))
	axis(2,labels=seq(min(lims[2,]),max(lims[2,]),(max(lims[2,])-min(lims[2,]))/4),at=seq(0,1,0.25))
	# contour(Mat,add=T,col="red",lwd=1,levels=0.5)

	## differentiate areas
	for(k in 1:length(MatList))
	{
		Mat_ = MatList[[k]]
		X = NULL
		Y = NULL
		Z = NULL
		for(i in seq(1,ncol(Mat_),10))
		{
			for(j in seq(1,ncol(Mat_),10))
			{
				z = Mat_[i,j]
				x = i
				y = j
				X = c(X,x)
				Y = c(Y,y)
				Z = c(Z,z)
			}
		}

		## add points
		X_ = (X-1)/(max(X)-1)
		Y_ = (Y-1)/(max(Y)-1)
		Z_ = (Z>0) + 2*(k-1)
		points(X_,Y_,pch=Z_,cex=Z)
	}
}

#
###

#######################
## UTILITY FUNCTIONS ##
#######################

## goal: functions for visualisation (figure 3 in the paper)

## .slice 
## goal: slice a vector x into a list of vectors of length specificied by y
.slice = function(x,y)
{
	xList = list()
	for(i in 1:length(y))
	{
		xList[[i]] = x[y[i-1] + 1:y[i]]
	}
	return(xList)
}

## .getmat ##
## goal: return matrix of function evaluated at grid points in variable limits 
## args: 
# lims - matrix   - limits of the two variables listed by rows
# func - function - function with two arguments
.getMat <- function(lims,func)
{
  
  ## compute 2D matrix
  res <- 100
  X <- apply(lims,1,function(X){seq(X[1],X[2],(X[2]-X[1])/res)}) # grid of pairwise combinations of x and y
  #
  ## fill output matrix - function value at each c(x,y)
  YMat <- matrix(0,ncol=res+1,nrow=res+1)
  for(i in 1:length(X[,1]))
  {
    for(j in 1:length(X[,2]))
    {
      YMat[i,j] <- func(c(X[i,1],X[j,2]))
    }
  }
  
  ## return matrix
  return(YMat)

}

## .heatmat ##
## goal: heatmap function for plotting 2D matrix
## args: 
# func - function - function with two arguments
# lims - vector   - min and max of the two variables
# labs - vector   - labels of the x and y axis
# main - vector   - main title of the plot
.heatMat <- function(YMat,lims=rbind(c(0,1),c(0,1)),labs=c("",""),main=c(""),axes=T)
{
  
  ## plot matrix
  #
  ## relative min and max to be matched to min and max color levels
  maxAbsMinMax <- max(abs(c(min(YMat),max(YMat))))
  #
  ## compute levels for color breaks
  levels <- seq(-maxAbsMinMax,maxAbsMinMax,2*maxAbsMinMax/1000)
  colorLevels <- rev(rainbow(1000,start=0,end=1,alpha=0.5))
  #
  ## heatmap
  image(YMat,breaks=levels,col=colorLevels,xaxt="n",yaxt="n",xlab=labs[1],ylab=labs[2],main=main)
  #
  ## add contour lines
  contour(YMat,add=T,col="black")
  #
  ## axes
  if(axes==T)
  {
  	for(i in 1:2){axis(side=i,at=c(0,.25,.5,.75,1),labels=round(c(0,.25,.5,.75,1)*diff(lims[i,])+lims[i,1],2))}
  }
  
}

## .plotMat ##
## goal: heatmap function for functions of two variables
## args: 
# lims - matrix   - limits of the two variables listed by rows
# func - function - function with two arguments
# labs - vector   - labels of the x and y axis
# main - vector   - main title of the plot
.plotMat <- function(lims,func,labs=c("",""),main=c(""),axes=T)
{
  YMat = .getMat(lims,func)
  .heatMat(YMat,lims,labs,main,axes)
}

#
###
