################
## analysis.R ##
################

## goal: analyse the results of the NODE fit to the simulated time series
## (i) visualise the fit to the time series
## (ii) compute and visualise direct effects on p.c. of hare and lynx (Jacobian)
## (iii) compute and visualise direct contributions on p.c. of hare and lynx (Geber method)
## (iv) visualise p.c. growth rate of hare and lynx (SLP bifurcation on hare and lynx density)

##############
## INITIATE ##
##############

## initiate
source("2-initiate.R")

## load ensemble
ensemble <- read.table("out/ensemble.csv",sep=";",header=T)

## filter only 25% best ensemble elements
threshold <- quantile(ensemble[,1],probs=0.25)
Omega_ <- ensemble[ensemble[,1]<threshold,-1]
print(dim(ensemble))
print(dim(Omega_))

#
###

######################################
## DIRECT EFFECTS AND CONTRIBUTIONS ##
######################################

## goal: compute and visualise (i) the fit of the NODE systems, (ii) direct effects and (iii) contributions to the p.c. growth rates of hare and lynx (SLPs)

## time steps and limits of the hare and lynx density
t <- seq(0,max(Y[,1]),0.1)
lims <- c(min(Y[,-1]),max(Y[,-1]))
#
## compute predicted states (NODE fit)
Ybar_ <- apply(Omega_,1,function(x)t(Ybar(t,as.matrix(Y[1,-1]),x))) # simulate NODE for all elements of the ensemble
Ybar_mean <- matrix(apply(Ybar_,1,mean),byrow=T,ncol=3) # ensemble average
Ybar_lo <- matrix(apply(Ybar_,1,function(x)quantile(probs=0.05,x)),byrow=T,ncol=3) # 5% quantile
Ybar_25 <- matrix(apply(Ybar_,1,function(x)quantile(probs=0.25,x)),byrow=T,ncol=3) # 25% quantile
Ybar_75 <- matrix(apply(Ybar_,1,function(x)quantile(probs=0.75,x)),byrow=T,ncol=3) # 75% quantile
Ybar_hi <- matrix(apply(Ybar_,1,function(x)quantile(probs=0.95,x)),byrow=T,ncol=3) # 95% quantile
#
## compute direct effects on p.c. growth rates (Jacobian)
Jbar_ <- apply(Omega_,1,function(x)Jbar(t,as.matrix(Y[1,-1]),x)) # compute Jacobian for all time steps and ensemble elements
Jbar_mean <- matrix(apply(Jbar_,1,mean),byrow=T,ncol=4) # ensemble average
Jbar_lo <- matrix(apply(Jbar_,1,function(x)quantile(probs=0.05,x)),byrow=T,ncol=4) # 5% quantile
Jbar_25 <- matrix(apply(Jbar_,1,function(x)quantile(probs=0.25,x)),byrow=T,ncol=4) # 25% quantile
Jbar_75 <- matrix(apply(Jbar_,1,function(x)quantile(probs=0.75,x)),byrow=T,ncol=4) # 75% quantile
Jbar_hi <- matrix(apply(Jbar_,1,function(x)quantile(probs=0.95,x)),byrow=T,ncol=4) # 95% quantile
#
## compute contributions to p.c. growth rates (Geber method)
Cbar_ <- apply(Omega_,1,function(x)Cbar(t,as.matrix(Y[1,-1]),x)) # compute contributions for all time steps and ensemble elements
Cbar_mean <- matrix(apply(Cbar_,1,mean),byrow=T,ncol=4) # ensemble average
Cbar_lo <- matrix(apply(Cbar_,1,function(x)quantile(probs=0.05,x)),byrow=T,ncol=4) # 5% quantile
Cbar_25 <- matrix(apply(Cbar_,1,function(x)quantile(probs=0.25,x)),byrow=T,ncol=4) # 25% quantile
Cbar_75 <- matrix(apply(Cbar_,1,function(x)quantile(probs=0.75,x)),byrow=T,ncol=4) # 75% quantile
Cbar_hi <- matrix(apply(Cbar_,1,function(x)quantile(probs=0.95,x)),byrow=T,ncol=4) # 95% quantile
#
## True effects and contributions (from LV model)
Ybar_true_ <- Ybar_true(t,c(1,1),c(1,.5,.5,1))
Jbar_true_ <- t(Jbar_true(t,c(1,1),c(1,.5,.5,1)))
Cbar_true_ <- t(Cbar_true(t,c(1,1),c(1,.5,.5,1)))

## plots
pdf("out/plots/figure_2.pdf")
#
par(mfrow=c(3,2))
#
## hare dynamics
i <- 2
plot(Y[,1],Y[,i],ylim=c(min(Ybar_lo[,i]),max(Ybar_hi[,i])),xlab="Time",ylab="Density",main="a.")
polygon(c(t,rev(t)),c(Ybar_lo[,i],rev(Ybar_hi[,i])),col=grey(0.5,0.2),border=NA) # 5%-95% quantile range
polygon(c(t,rev(t)),c(Ybar_25[,i],rev(Ybar_75[,i])),col=grey(0.5,0.4),border=NA) # 25-75% quantile range
lines(t,Ybar_mean[,i],col="red") # NODE predictions
#
## lynx dynamics
i <- 3
plot(Y[,1],Y[,i],ylim=c(min(Ybar_lo[,i]),max(Ybar_hi[,i])),xlab="Time",ylab="Density",main="b.")
polygon(c(t,rev(t)),c(Ybar_lo[,i],rev(Ybar_hi[,i])),col=grey(0.5,0.2),border=NA) # 5%-95% quantile range
polygon(c(t,rev(t)),c(Ybar_25[,i],rev(Ybar_75[,i])),col=grey(0.5,0.4),border=NA) # 25-75% quantile range
lines(t,Ybar_mean[,i],col="red") # NODE predictions
#
## direct effects on hare p.c. growth
i <- 1
plot(t,Jbar_mean[,i],ylim=c(min(Jbar_lo[,c(i,i+1)]),max(Jbar_hi[,c(i,i+1)])),cex=0,xlab="Time",ylab="Effect",main="c.")
#
polygon(c(t,rev(t)),c(Jbar_lo[,i],rev(Jbar_hi[,i])),col=grey(0.5,0.2),border=NA) # 5%-95% quantile range
polygon(c(t,rev(t)),c(Jbar_25[,i],rev(Jbar_75[,i])),col=grey(0.5,0.4),border=NA) # 25-75% quantile range
lines(t,Jbar_mean[,i],col="black") # intra-specific effect
lines(t,Jbar_true_[,i],lty=2,col="black") # intra-specific effect (ground truth)
#
i <- 2
polygon(c(t,rev(t)),c(Jbar_lo[,i],rev(Jbar_hi[,i])),col=grey(0.5,0.2),border=NA) # 5%-95% quantile range
polygon(c(t,rev(t)),c(Jbar_25[,i],rev(Jbar_75[,i])),col=grey(0.5,0.4),border=NA) # 25-75% quantile range
lines(t,Jbar_mean[,i],col="red") # inter-specific effect
lines(t,Jbar_true_[,i],lty=2,col="red") # inter-specific effect (ground truth)
#
## direct effects on lynx p.c. growth
i <- 3
plot(t,Jbar_mean[,i],ylim=c(min(Jbar_lo[,c(i,i+1)]),max(Jbar_hi[,c(i,i+1)])),cex=0,xlab="Time",ylab="Effect",main="d.")
#
polygon(c(t,rev(t)),c(Jbar_lo[,i],rev(Jbar_hi[,i])),col=grey(0.5,0.2),border=NA) # 5%-95% quantile range
polygon(c(t,rev(t)),c(Jbar_25[,i],rev(Jbar_75[,i])),col=grey(0.5,0.4),border=NA) # 25-75% quantile range
lines(t,Jbar_mean[,i],col="red") # inter-specific effect
lines(t,Jbar_true_[,i],lty=2,col="red") # inter-specific effect (ground truth)
#
i <- 4
polygon(c(t,rev(t)),c(Jbar_lo[,i],rev(Jbar_hi[,i])),col=grey(0.5,0.2),border=NA) # 5%-95% quantile range
polygon(c(t,rev(t)),c(Jbar_25[,i],rev(Jbar_75[,i])),col=grey(0.5,0.4),border=NA) # 25-75% quantile range
lines(t,Jbar_mean[,i],col="black") # intra-specific effect
lines(t,Jbar_true_[,i],lty=2,col="black") # intra-specific effect (ground truth)
#
## contributions to hare p.c. growth
i <- 1
plot(t,Cbar_mean[,i],ylim=c(min(Cbar_lo[,c(i,i+1)]),max(Cbar_hi[,c(i,i+1)])),cex=0,xlab="Time",ylab="Contribution",,main="e.")
#
polygon(c(t,rev(t)),c(Cbar_lo[,i],rev(Cbar_hi[,i])),col=grey(0.5,0.2),border=NA) # 5%-95% quantile range
polygon(c(t,rev(t)),c(Jbar_25[,i],rev(Jbar_75[,i])),col=grey(0.5,0.4),border=NA) # 25-75% quantile range
lines(t,Cbar_mean[,i],col="black") # intra-specific contribution
lines(t,Cbar_true_[,i],lty=2,col="black") # intra-specific contribution (ground truth)
#
i <- 2
polygon(c(t,rev(t)),c(Cbar_lo[,i],rev(Cbar_hi[,i])),col=grey(0.5,0.2),border=NA) # 5%-95% quantile range
polygon(c(t,rev(t)),c(Cbar_25[,i],rev(Cbar_75[,i])),col=grey(0.5,0.4),border=NA) # 25-75% quantile range
lines(t,Cbar_mean[,i],col="red") # inter-specific contribution
lines(t,Cbar_true_[,i],lty=2,col="red") # inter-specific contribution (ground truth)
#
## contributions to lynx p.c. growth
i <- 3
plot(t,Cbar_mean[,i],ylim=c(min(Cbar_lo[,c(i,i+1)]),max(Cbar_hi[,c(i,i+1)])),cex=0,xlab="Time",ylab="Contribution",,main="f.")
#
polygon(c(t,rev(t)),c(Cbar_lo[,i],rev(Cbar_hi[,i])),col=grey(0.5,0.2),border=NA)
polygon(c(t,rev(t)),c(Cbar_25[,i],rev(Cbar_75[,i])),col=grey(0.5,0.4),border=NA) # 25-75% quantile range
lines(t,Cbar_mean[,i],col="red") # inter-specific contribution
lines(t,Cbar_true_[,i],lty=2,col="red") # inter-specific contribution (ground truth)
#
i <- 4
polygon(c(t,rev(t)),c(Cbar_lo[,i],rev(Cbar_hi[,i])),col=grey(0.5,0.2),border=NA)
polygon(c(t,rev(t)),c(Cbar_25[,i],rev(Cbar_75[,i])),col=grey(0.5,0.4),border=NA) # 25-75% quantile range
lines(t,Cbar_mean[,i],col="black") # intra-specific contribution
lines(t,Cbar_true_[,i],lty=2,col="black") # intra-specific contribution (ground truth)
#
par(mfrow=c(1,1))
#
dev.off()

#
###

#################################
## VISUALISE P.C. GROWTH RATES ##
#################################

## goal: (iv) visualise SLP approximations of the p.c. growth rates of hare and lynx

## evaluate SLP over all combinations of hare and lynx density
SLPMat_ <- apply(Omega_,1,function(X)SLPMat(lims,lims,X)) # SLP evaluated at every pairwise combination of hare in lynx in lims
SLPMat_mean <- matrix(apply(SLPMat_,1,mean),byrow=T,ncol=2) # ensemble average

## true per-capita growth rates
r1_true <- function(state){1-0.5*state[2]}
r2_true <- function(state){0.5*state[1]-1}

## plots
pdf("out/plots/figure_3.pdf")
#
par(mfrow=c(2,2))
#
## scale data (to add it to heatmap below)
Y_ <- apply(Y[,-1],2,function(x)(x-lims[1])/(lims[2]-lims[1]))
#
## hare p.c. growth rate
.heatmat(t(matrix(SLPMat_mean[,1],ncol=101)),lims=rbind(lims,lims),labs=c("Hare density","Lynx density"),main="a.")
points(Y_[,1],Y_[,2],pch=16) # data points
#
## hare p.c. growth rate (ground truth)
.heatmap(func=r1_true,lims=rbind(lims,lims),labs=c("Hare density","Lynx density"),main="b.")
points(Y_[,1],Y_[,2],pch=16) # data points
#
## lynx p.c. growth rate
.heatmat(t(matrix(SLPMat_mean[,2],ncol=101)),lims=rbind(lims,lims),labs=c("Hare density","Lynx density"),main="c.")
points(Y_[,1],Y_[,2],pch=16) # data points
#
## lynx p.c. growth rate (ground truth)
.heatmap(func=r2_true,lims=rbind(lims,lims),labs=c("Hare density","Lynx density"),main="d.")
points(Y_[,1],Y_[,2],pch=16) # data points
#
par(mfrow=c(1,1))
#
dev.off()

#
###
