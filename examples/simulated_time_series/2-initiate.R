################
## initiate.R ##
################

## goal: initiate the NODE system and the simulated time series

##############
## INITIATE ##
##############

## load ODE library
library(deSolve)

## source functions
source("1-functions.R")

#
###

###########################
## SIMULATED TIME SERIES ##
###########################

## goal: generate the simulated time series from the Lotka-Volterra equations

## initiate LV model
Y_0 <- c(1,1) # initial state
param <- c(1,.5,.5,1) # parameters
#
## simulate LV ODE system
Y <- data.frame(Ybar_true(seq(0,7,.5),Y_0,param)) 
Ymax <- max(Y[,-1]); # to scale the system in the SLPs

#
###
