#################
## functions.R ##
#################

## goal: declare all functions for the NODE analysis of the real time series

#################
## NODE SYSTEM ##
#################

## goal: define the NODE system

## ANN dimensions and indices
D <- 10 # number of hidden neurones
d0 <- 1:D # intermediate-to-output weights indices
d1 <- 1:D + D # time weight indices
d2 <- 1:D + 2*D # hare density weight indices
d3 <- 1:D + 3*D # lynx density weight indices
d4 <- 1:D + 4*D # intermediate nodes bias indices
d5 <- 5*D + 1 # output node bias index
d <- d5 # total newtork size

## SLP ##
## goal: evaluate the SLP given the input vector X and weight vector Omega
## args:
## @ X - vector - vector of input variables
## @ Omega - vector - vector of biases and weights for the SLP
SLP <- function(X,Omega)
{
  return(Omega[d0]%*%f_sigma(X[1]*Omega[d1]/30 + X[2]*Omega[d2] + X[3]*Omega[d3] + Omega[d4])+Omega[d5])
}

## SLPMat ##
## goal: evaluate SLPs over a matrix of pairwise combinations of inputs
## args:
## @ xlims, ylims - vector - vector of min and max values between which to evaluate the SLP
## @ param - vector - vector of biases and weights of the SLPs
SLPMat <- function(t,xlims,ylims,param)
{
  ##
  x_ <- seq(xlims[1],xlims[2],(xlims[2]-xlims[1])/100)
  y_ <- seq(ylims[1],ylims[2],(ylims[2]-ylims[1])/100)
  
  ##
  return(apply(pair(x_,y_),1,function(X)c(SLP(c(t,X),param[1:d]),SLP(c(t,X),param[1:d+d]))))
}

## f_sigma - sigmoid ##
## goal: SLP sigmoid activation function
## args:
## @ x - vector - vector of inputs
f_sigma <- function(x)
{
  return(1/(1+exp(-x)))
}

## dSLPdx - sigmoid ##
## goal: gradient of the SLP with respect to each input (valid for the sigmoid activation function only)
# args:
## @ X - vector - vector of inputs
## @ Omega - vector - vector of biases and weights of the SLP
dSLPdx <- function(X,Omega)
{

  ## weighted sum of inputs 
  u <- X[1]*Omega[d1]/30 + X[2]*Omega[d2] + X[3]*Omega[d3] + Omega[d4]

  ## gradient wtr inputs (hare and lynx density only)
  return(c(
    (Omega[d0]*Omega[d2])%*%((1-f_sigma(u))*f_sigma(u)),
    (Omega[d0]*Omega[d3])%*%((1-f_sigma(u))*f_sigma(u))
  ))

}

## Ydot ##
## goal: ordinary differential equation part of the NODE system
## args: 
## @ t - float - time step at which to evaluate the differential equations (not useful here)
## @ state - vector - vector of states (hare and lynx density)
## @ param - vector - vector of biases and weights of the SLPs
Ydot <- function(t,state,param)
{
  return(list(c(SLP(c(t,state),param[1:d])*state[1], SLP(c(t,state),param[1:d+d])*state[2])))
}

## Ybar ##
## goal: predictive function to simulate the NODE system
## args: 
## @ t - float - time step at which to evaluate the differential equations
## @ state_0 - vector - vector of initial states
## @ param - vector - vector of biases and weights of the SLPs
Ybar <- function(t,state_0,param)
{
  return(ode(y=state_0, times=t, func=Ydot, parms=param, method="ode45", maxsteps=100))
}

## error ##
## goal: compute the error in NODE predictions (Ordinary Least Squares)
## args:
# param - vector - vector of biases and weights of the SLPs
error <- function(param)
{
  return(1/(nrow(Y)*2)*sum((Y[,-1]-Ybar(Y[,1],as.matrix(Y[1,-1]),param)[,-1])^2))
}

## J ##
## goal: compute the Jacobian of both SLPs at a given time step
## args:
## @ t - float - time step at which to compute the Jacobian
## @ state - vector - vector of states
## @ param - vector - vector of biases and weights of the SLPs
J <- function(t,state,param)
{
  return(c(dSLPdx(c(t,state),param[1:d]),dSLPdx(c(t,state),param[1:d+d])))
}

## Jbar ##
## goal: compute the Jacobian of each SLPs across several time steps
## args:
## @ t - vector - vector of time steps at which to compute the Jacobian
## @ state_0 - vector - vector of initial states
## @ param - vector - vector of biases and weights of the SLPs
Jbar <- function(t,state_0,param)
{
  return(apply(Ybar(t,state_0,param),1,function(X)J(X[1],X[-1],param)))
}

## C ##
## goal: compute the contribution of each variable to both SLPs at a given time step (chaing rule/Geber method)
## args:
## @ t - float - time step at which to compute the Jacobian
## @ state - vector - vector of states
## @ param - vector - vector of biases and weights of the SLPs
C <- function(t,state,param)
{
  return(rep(unlist(Ydot(t,state,param)),2)*J(t,state,param))
}

## Cbar ##
## goal: compute the contiburions of each variable to SLPs across several time steps (chain rule/Geber method)
## args:
## @ t - vector - vector of time steps at which to compute the contributions
## @  state_0 - vector - vector of initial states
## @  param - vector - vector of biases and weights of the SLPs
Cbar <- function(t,state_0,param)
{
  return(apply(Ybar(t,state_0,param),1,function(X)C(X[1],X[-1],param)))
}

#
###

########################
## PLOTTING FUNCTIONS ##
########################

## goal: functions for visualisation (figure 3 in the paper)

## .heatmat ##
## goal: heatmap function for plotting 2D matrix
## args: 
## @ func - function - function with two arguments
## @ lims - vector - min and max of the two variables
## @ labs - vector - labels of the x and y axis
## @ main - vector - main title of the plot
.heatmat <- function(YMat,lims=rbind(c(0,1),c(0,1)),labs=c("",""),main=c(""))
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
  for(i in 1:2){axis(side=i,at=c(0,.25,.5,.75,1),labels=round(c(0,.25,.5,.75,1)*diff(lims[i,])+lims[i,1],2))}
  
}

## .heatmap ##
## goal: heatmap function for functions of two variables
## args: 
## @ lims - matrix - limits of the two variables listed by rows
## @ func - function - function with two arguments
## @ labs - vector - labels of the x and y axis
## @ main - vector - main title of the plot
.heatmap <- function(lims,func,labs=c("",""),main=c(""))
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
  
  ## plot matrix
  .heatmat(YMat,lims,labs,main)
}

#
###

#######################
## UTILITY FUNCTIONS ##
#######################

## goal: utility functions 

## pair ##
## goal: create a 2 column matrix of pairwise combinations of all elements of two vectors
## args:
## @ x, y - vector - vectors to combine pairwise
pair <- function(x,y)
{
  return(cbind(as.vector(matrix(x,ncol=length(x),nrow=length(y),byrow=T)),y))
}

#
###
