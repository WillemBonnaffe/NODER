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
# 07-01-2021 - created version 0.4
# 07-01-2021 - replaced polynomial process model by ANN process model 
# 08-01-2021 - created version 0.5
# 08-01-2021 - decoupled the interpolation from the fitting of the process model

## to do
# - compute r2 for each time series
# - visualisation of the interpolation 

##############
## INITIATE ##
##############

## goal: initiate the RMT algorithm

# ## libraries
# library(Rcpp)
# sourceCpp("cpp/DEMCpp_v0.1.cpp")
# sourceCpp("cpp/DEMCO_v0.5.cpp")

## TS
TS = read.table("TS.csv",sep=";",header=T)

## focus on fortis for now
TS = TS[,c("year","rain","small.seeds","large.seeds","fortis.N","fortis.PC1.beak","fortis.PC2.beak")]

## transform real variables
TS[,c("rain","small.seeds","large.seeds","fortis.N")] = log(TS[,c("rain","small.seeds","large.seeds","fortis.N")])

#
###

###############
## FUNCTIONS ##
###############

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

logPost_o = function(t,Y,Omega_o)
{
	## predict
	Ybar_o = f_o.eval(t,Omega_o)

	## residuals
	res_o = Y - Ybar_o

	## Likelihood
	# logLik_o   = sum(log(dnorm(res_o,0,sd1_o)))
	logLik_o   = sum(-0.5*log(2*pi*(sd1_o^2)) - 0.5*(res_o^2)/(sd1_o^2))

	## prior
	# logPrior_o = sum(log(dnorm(Omega_o,0,sd2_o)))
	logPrior_o = sum(-0.5*log(2*pi*(sd2_o^2)) - 0.5*(Omega_o^2)/(sd2_o^2))

	## posterior
	logPost_o  = logLik_o + logPrior_o

	## terminate
	return(logPost_o)
}

ddOmega_o.logPost_o = function(t,Y,Omega_o)
{
	## predict
	Ybar_o = f_o.eval(t,Omega_o)

	## derivative of residuals
	res_o = (Y - Ybar_o)
	ddOmega_o.res_o = -ddOmega_o.f_o.eval(t,Omega_o)

	## derivative of likelihood
	ddOmega_o.logLik_o = -ddOmega_o.res_o%*%res_o/(sd1_o^2)

	## derivative of prior
	ddOmega_o.logPrior_o = -Omega_o/(sd2_o^2)

	## derivative of posterior
	ddOmega_o.logPost_o = ddOmega_o.logLik_o + ddOmega_o.logPrior_o

	## terminate
	return(ddOmega_o.logPost_o)
}

argmin.logPost_o = function(t,Y,Omega_o)
{
	error_     = function(x) -logPost_o(t,Y,x)
	graderror_ = function(x) -ddOmega_o.logPost_o(t,Y,x)
	Omega_o  = optim(par=Omega_o,
					 fn=error_,
					 gr=graderror_,
					 method="BFGS"#,
					 # control=list("trace"=1,"REPORT"=1,"maxit"=100)
					 )$par
	return(Omega_o)
}

## functions for the process model
f_sigma     = function(x) 1/(1+exp(-x))
ddu.f_sigma = function(x) f_sigma(x) * (1 - f_sigma(x))
# f_sigma     = function(x) x
# ddu.f_sigma = function(x) 1
ddt.f_p = function(x,Omega) 
{	
	Omega = matrix(Omega,nrow=W_p)
	return(t(Omega[,1])%*%f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x))
}
ddOmega_p.ddt.f_p = function(x,Omega)
{
	Omega = matrix(Omega,nrow=W_p)
	ddOmega_p1. = f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x)
	ddOmega_p2. = Omega[,1]*1*ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x)
	ddOmega_p3. = Omega[,1]*rep(1,W_p)%*%t(x)*as.vector(ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x))
	return(c(ddOmega_p1.,ddOmega_p2.,ddOmega_p3.))
}
ddx.ddt.f_p = function(x,Omega)
{
	Omega = matrix(Omega,nrow=W_p)
	ddx. = Omega[,1]%*%(Omega[,-c(1:2)]*as.vector(ddu.f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%x)))
	return(ddx.)
}
ddt.f_p.eval           = function(X,Omega) apply(X,1,function(x) ddt.f_p(x,Omega))
ddx.ddt.f_p.eval       = function(X,Omega) apply(X,1,function(x) ddx.ddt.f_p(x,Omega))
ddOmega_p.ddt.f_p.eval = function(X,Omega) apply(X,1,function(x) ddOmega_p.ddt.f_p(x,Omega))

logPost_p = function(Ybar_o,ddt.Ybar_o,Omega_p)
{
	## predict
	ddt.Ybar_p = ddt.f_p.eval(Ybar_o,Omega_p) 

	## residuals
	res_p = ddt.Ybar_o - ddt.Ybar_p

	## Likelihood
	logLik_p   = sum(-0.5*log(2*pi*(sd1_p^2)) - 0.5*(res_p^2)/(sd1_p^2))

	## prior
	# logPrior_o = sum(log(dnorm(Omega_o,0,sd2_o)))
	logPrior_p = sum(-0.5*log(2*pi*(sd2_p^2)) - 0.5*(Omega_p^2)/(sd2_p^2))

	## posterior
	logPost_p  = logLik_p + logPrior_p

	## terminate
	return(logPost_p)
}

ddOmega_p.logPost_p = function(Ybar_o,ddt.Ybar_o,Omega_p)
{
	## predict
	ddt.Ybar_p = ddt.f_p.eval(Ybar_o,Omega_p) 

	## derivative of residuals
	res_p = ddt.Ybar_o - ddt.Ybar_p
	ddOmega_p.res_p = -ddOmega_p.ddt.f_p.eval(Ybar_o,Omega_p)

	## derivative of likelihood
	ddOmega_p.logLik_p = -ddOmega_p.res_p%*%res_p/(sd1_p^2)

	## derivative of prior
	ddOmega_p.logPrior_p = -Omega_p/(sd2_p^2)

	## derivative of posterior
	ddOmega_p.logPost_p = ddOmega_p.logLik_p + ddOmega_p.logPrior_p

	## terminate
	return(ddOmega_p.logPost_p)
}

argmin.logPost_p = function(Ybar_o,ddt.Ybar_o,Omega_p)
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

##
resTab  = NULL
for( lambda in seq(0.001,0.03,0.001))
{

## initiate model
n     = nrow(TS)
n_    = n*2
N     = ncol(TS) - 1
W_o   = 10
W_p   = 10
sd1_o = 0.01
sd2_o = 0.1
sd1_p = 0.01
sd2_p = lambda # 0.015 
alpha = 10
t     = TS[,1]
Y     = TS[,-1] 
t     = (t-min(t))/(max(t)-min(t))*alpha
mean_ = apply(Y,2,mean)
sd_   = apply(Y,2,sd)
Y     = t(apply(Y,1,function(x)(x-mean_)/sd_)) # standardise Y
t_    = seq(min(t),max(t),(max(t)-min(t))/n_)

## fit observation model
# load("Omega_oList.RData")
# Omega_oList = list()
# for(i in 1:N)
# {
# 	## BFGS
# 	Omega_o = rnorm(3*W_o,0,sd2_o)
# 	Omega_o = argmin.logPost_o(t,Y[,i],Omega_o) 
# 	Omega_oList[[i]] = Omega_o
# }
# save(Omega_oList,file="Omega_oList.RData")

## get time series
Ybar_o = NULL
ddt.Ybar_o = NULL
for(i in 1:N)
{
	Omega_o     = Omega_oList[[i]]
	Ybar_o      = cbind(Ybar_o,f_o.eval(t_,Omega_o))
	ddt.Ybar_o  = cbind(ddt.Ybar_o,ddt.f_o.eval(t_,Omega_o))
}

## fit process model
Omega_pList = list()
for(i in 1:N)
{
	## BFGS
	Omega_p = rnorm((2+N)*W_p,0,sd2_p)
	Omega_p = argmin.logPost_p(Ybar_o,ddt.Ybar_o[,i],Omega_p) 
	Omega_pList[[i]] = Omega_p
}

## get process model
ddt.Ybar_p = NULL
ddx.ddt.Ybar_pList = list()
for(i in 1:N)
{
	Omega_p = Omega_pList[[i]]
	ddt.Ybar_p  = cbind(ddt.Ybar_p,ddt.f_p.eval(Ybar_o,Omega_p))
	ddx.ddt.Ybar_pList[[i]]  = t(ddx.ddt.f_p.eval(Ybar_o,Omega_p))
}

## get process residuals
res_p = ddt.Ybar_o - ddt.Ybar_p 
sd.res_p      = apply(res_p,2,sd)
sd.ddt.Ybar_o = apply(ddt.Ybar_o,2,sd)
r2            = 1 - sd.res_p^2/sd.ddt.Ybar_o^2
resTab = rbind(resTab,c(lambda,r2))

## DEMC
# Omega_o  = rnorm(N*3*W_o,0,sd2_o)
# error_     = function(x) logPost_o(t,Y[,1],x)
# chainList = DEMCOpp(list("dTarget"=error_,
# 						 "Theta_0"=Omega_o,
# 						 "epsilon"=0.001,
# 						 "lambda"=100,
# 						 "nIt"=1000))
# Omega_o = chainList$chain[which.max(chainList$chain[,1]),-1]

## DEMC
# Omega_o  = rnorm(N*3*W_o,0,sd2_o)
# error_     = function(x) logPost_o(t,Y[,1],x)
# chainList = DEMCpp(list("dTarget"=error_,"Theta_0"=Omega_o,"epsilon"=0.001,"nIt"=10000))
# Omega_o = chainList$chain[which.max(chainList$chain[,1]),-1]

}

#
###


##############
## ANALYSIS ##
##############

pdf("results.pdf")

## visualise change in residuals
plot(resTab[,1],resTab[,2],ylim=c(0,1),cex=0)
for(i in 1:N) lines(resTab[,1],resTab[,i+1],col=rainbow(N)[i])

## visualise time series
par(mfrow=c(N,N),mar=c(1,1,1,1),oma=c(1,1,1,1))
#
for(i in 1:N)
{
	for(j in 1:N)
	{
		if(i == j)
		{
			plot(-1:1,-1:1,cex=0,xaxt="n",yaxt="n")
			text(0,0,colnames(Y)[i])
		} else
		{
			plot(Y[,i],Y[,j],xaxt="n",yaxt="n",type="l")
			lines(Ybar_o[,i],Ybar_o[,j],col="red")
		}
	}
}
#
par(mfrow=c(1,1),mar=c(5,4,4,2)+.1,oma=c(1,1,1,1)*0)


## visualise fit
layout(mat=cbind(c(1,2,3,4),c(5,6,7,8)))
par(mar=c(4,4,1,1)+0.1,oma=c(1,1,1,1))
colVect = rainbow(N)
for(i in 1:N)
{
	plot(t,Y[,i],type="b",ylim=c(-3,3))
	lines(t_,Ybar_o[,i],col=rainbow(N)[i])
	lines(c(min(t),max(t)),c(0,0),lty=2)
	#
	plot(t_,ddt.Ybar_o[,i],col="blue",ylim=c(-4,4),type="l")
	lines(t_,ddt.Ybar_p[,i],col="green")
	lines(c(min(t),max(t)),c(0,0),lty=2)
	#
	plot(t_,ddx.ddt.Ybar_pList[[i]][,1],ylim=c(-3,3),type="l")
	for(j in 1:N)
	{
		lines(t_,ddx.ddt.Ybar_pList[[i]][,j],col=rainbow(N)[j])
	}
	lines(c(min(t),max(t)),c(0,0),lty=2)
	#
	plot(t_,ddt.Ybar_o[,1]*ddx.ddt.Ybar_pList[[i]][,1],ylim=c(-3,3),type="l")
	for(j in 1:N)
	{
		lines(t_,ddt.Ybar_o[,j]*ddx.ddt.Ybar_pList[[i]][,j],col=rainbow(N)[j])
	}
	lines(c(min(t),max(t)),c(0,0),lty=2)

}
par(mfrow=c(1,1),mar=c(5,4,4,2)+.1,oma=c(1,1,1,1)*0)

## compute contribution matrix
EMat = matrix(0,nrow=N,ncol=N)
for(i in 2:N) # skip rain
{
	for(j in 1:N)
	{
		EMat[i,j] = mean(ddx.ddt.Ybar_pList[[i]][,j])
	}
	# EMat[i,] = EMat[i,]/sum(EMat[i,])
}

## compute contribution matrix
CMat = matrix(0,nrow=N,ncol=N)
for(i in 2:N) # skip rain
{
	for(j in 1:N)
	{
		CMat[i,j] = sum((ddt.Ybar_o[,j]*ddx.ddt.Ybar_pList[[i]][,j])^2)
	}
	CMat[i,] = CMat[i,]/sum(CMat[i,])
}

## visualise contribution matrix
# image(CMat,col=rev(rainbow(100,start=0.2,end=0.8)))
.plot.DIN(EMat,CMat*(CMat>0.1),colnames(Y))

# ## plot contribution matrix
# layout(mat=cbind(c(1,2,3),c(5,6,7)))
# colVect = rainbow(N)
# for(i in 1:N)
# {
# 	#
# 	plot(t_,ddt.Ybar_o[,1]*ddx.ddt.Ybar_pList[[i]][,1],ylim=c(-3,3),type="l")
# 	for(j in 1:N)
# 	{
# 		lines(t_,ddt.Ybar_o[,j]*ddx.ddt.Ybar_pList[[i]][,j],col=rainbow(N)[j])
# 	}
# 	lines(c(min(t),max(t)),c(0,0),lty=2)
# }
# par(mfrow=c(1,1))

## visualise phase space

## residuals checks
par(mfrow=c(2,2))
for(i in 1:N)
{
	# res_o = Y - Ybar_o[s,]
	# barplot(res_o[,i])
	barplot(res_p[,i])
	legend("topright",legend=round(r2[i],2))
}
par(mfrow=c(1,1))

dev.off()

#
###
