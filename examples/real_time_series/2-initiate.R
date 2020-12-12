################
## initiate.R ##
################

## goal: initiate the NODE system and the real time series

##############
## INITIATE ##
##############

## load ODE library
library(deSolve)

## source functions
source("1-functions.R")

#
###

#################
## TIME SERIES ##
#################

## goal: load the time series from hare lynx pelt counts from the Hudson Bay company records (1904-1934)

## read data
Y <- read.table(file = "hare_lynx_data.csv", sep=";", header=T) 
Y <- Y[60:90,] # last 31 years
Y[,1] <- Y[,1] - min(Y[,1]) # first year (1904) should be 0
Y[,-1] <- apply(Y[,-1],2,function(x)x/10) # change Y unit to hundreds of thousands

#
###
