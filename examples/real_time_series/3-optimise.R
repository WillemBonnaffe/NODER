################
## optimise.R ##
################

## goal: fit an ensemble of NODE systems to the real time series

##############
## INITIATE ##
##############

## initiate
source("2-initiate.R")

#
###

#################
## OPTIMSATION ##
#################

## goal: optimise the biases and weights of the SLP given the simulated time series following an ensemble approach

## initialise ensemble
ensemble <- NULL
nElements <- 200
nIt <- 150
#
## for each element of the ensemble
for(i in 1:nElements)
{
  ## display progress
  message(paste(i,"/",nElements))
  #
  ## initialise biases and weights
  Omega <- rnorm(d*2,0,0.001)
  #
  ## run BFGS optimisation
  optimised <- optim(par=Omega, fn=error, method="BFGS",control=list("trace"=1,"REPORT"=10,"maxit"=nIt))
  #
  ## add element to ensemble
  ensemble <- rbind(ensemble, c(optimised$value,optimised$par))
}

## store ensemble
write.table(ensemble,"out/ensemble.csv",sep=";",row.names=FALSE,col.names=FALSE,quote=FALSE,append=T)

#
###
