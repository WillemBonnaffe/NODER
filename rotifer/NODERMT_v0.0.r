##########################
## NODE_RMT_functions.r ##
##########################

## goal: perfrom residual minimisation training on the finch time series

## method:
# 1. interpolate the time series with sin ANN functions
# 2. estimate linear and non-linear coupling between the time series

## versions:
# 17-12-2020 - created version 0.0

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
n         = nrow(TS)
N         = ncol(TS) - 1
lambda_o  = 0.2
lambda_p  = 0.2
alpha     = 10
t         = TS[,1]
Y         = TS[,-1] 
## std t
t         = (t-min(t))/(max(t)-min(t))*alpha
## std Y
mean_     = apply(Y,2,mean)
sd_       = apply(Y,2,sd)
Y         = t(apply(Y,1,function(x)(x-mean_)/sd_)) # standardise Y

## check time series
plot(TS)


## functions for the observation model
f_o = function(x,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	return(t(Omega[,1])%*%sin(pi*(x*Omega[,2] + Omega[,3])))
}
f_o.eval = function(x,Omega)
{
	apply(t(x),2,function(x)f_o(x,Omega))
}
ddx.f_o = function(x,Omega)
{	
	Omega = matrix(Omega,ncol=3)
	return(t(Omega[,1])%*%(pi*Omega[,2]*cos(pi*(x*Omega[,2] + Omega[,3]))))
}
ddx.f_o.eval = function(x,Omega)
{
	apply(t(x),2,function(x)ddx.f_o(x,Omega))
}
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
	apply(t(x),2,function(x)ddOmega_o.f_o(x,Omega))
}
error.f_o = function(Y,Omega)
{
	mean((Y-f_o.eval(t,Omega))^2) + lambda_o*sum(Omega^2)
}
ddOmega_o.error.f_o = function(Y,Omega)
{
	as.vector(1/n*2*(-ddOmega_o.f_o.eval(t,Omega))%*%(Y-f_o.eval(t,Omega)) + lambda_p*2*Omega)
}

## functions for the process model
ddx.f_p             = function(X,Omega) X%*%Omega
ddOmega_p.ddx.f_p   = function(X,Omega) X
error.f_p           = function(ddx.f_o,Omega) mean((ddx.f_o - ddx.f_p(X_p,Omega))^2) + lambda_p*sum(Omega^2)
ddOmega_p.error.f_p = function(ddx.f_o,Omega) 1/n*2*t(-ddOmega_p.ddx.f_p(X_p,Omega))%*%(ddx.f_o - ddx.f_p(X_p,Omega)) + lambda_p*2*Omega

#
###

###############
## FIT MODEL ##
###############

## goal: fit the model to the time series

## fit observation model
Omega_o = NULL
Ybar_o  = NULL
dYbar_o = NULL
for(i in 1:N)
{
	## train
	Omega_o_ = rnorm(3*10)
	error_   = function(x)error.f_o(Y[,i],x)
	grad_    = function(x)ddOmega_o.error.f_o(Y[,i],x)
	Omega_o_ = optim(par=Omega_o_,fn=error_,gr=grad_,method="BFGS")$par

	## predict
	t_       = seq(min(t),max(t),(max(t)-min(t))/100)
	Ybar_o_  = f_o.eval(t_,Omega_o_)
	dYbar_o_ = ddx.f_o.eval(t_,Omega_o_)

	## store predictions
	Omega_o = cbind(Omega_o,Omega_o_)
	Ybar_o  = cbind(Ybar_o,Ybar_o_)
	dYbar_o = cbind(dYbar_o,dYbar_o_)
}

plot(data.frame(Ybar_o))

## fit process model
Omega_p = NULL
dYbar_p = NULL
# dYbar_p_num = NULL # for debugging
PIM = function(x) # pairwise interaction matrix function
{ 
	PIM_ = x%*%t(x)
	PIM_diag     = diag(PIM_)
	PIM_lowertri = as.vector(PIM_)[as.vector(lower.tri(PIM_))]
	return(c(PIM_diag,PIM_lowertri))
}
X_p     = cbind(1,Ybar_o[,1],Ybar_o[,2],Ybar_o[,3],Ybar_o[,4],Ybar_o[,5],Ybar_o[,6],t(apply(Ybar_o,1,PIM)))

for(i in 1:N) # skip rain
{
  ## train
	Omega_p_    = rnorm(ncol(X_p))
	error_      = function(x){error.f_p(dYbar_o[,i],x)}
	grad_       = function(x){ddOmega_p.error.f_p(dYbar_o[,i],x)}
	Omega_p_    = optim(par=Omega_p_,fn=error_,gr=grad_,method="BFGS")$par
	# Omega_p_num_ = optim(par=Omega_p_,fn=error_,method="BFGS")$par # for debugging

	## predict
	dYbar_p_    = ddx.f_p(X_p,Omega_p_)
	# dYbar_p_num_ = ddx.f_p(X_p,Omega_p_num_) # for debugging

	## store predicition
	Omega_p     = cbind(Omega_p,Omega_p_)
	dYbar_p     = cbind(dYbar_p,dYbar_p_)
	# dYbar_p_num = cbind(dYbar_p_num,dYbar_p_num_) # for debugging
}

#
###

##############
## ANALYSIS ##
##############

pdf("results.pdf")

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
